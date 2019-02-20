##Timothy R. Fallon 2019
##Golden Gate in-silico cloning automator.
# -*- coding: UTF-8 -*-

import tornado
import tornado.ioloop
import tornado.web
import os, uuid
import os.path
import io
import Bio
import Bio.SeqIO
import Bio.Seq
import Bio.SeqRecord
import Bio.Alphabet
import glob
import re
import logging
import base64
import zipfile
import sys
import auto_golden_gate_functions
 
TCP_PORT = 5686

class Userform(tornado.web.RequestHandler):
    def get(self):
        self.render("fileuploadform.html")
  
class Upload(tornado.web.RequestHandler):
    def post(self):
    ##Handle clicking of the "submit" button
        os.chdir(os.path.dirname(os.path.realpath(sys.argv[0]))) ##Change to the directory of the script.
   
        self.write("<html><font face=\"courier\">")
        self.write("<div style =\"width: 5000px; overflow:scroll\">") ##Some tricks to have scrolling text.
        self.write("parsing...")
        plasmid_name = self.get_argument('plasmid_name', '')
        text_box = self.get_argument('text_box','').strip().replace('\n','').replace(' ','')
        record_id = self.get_argument('record_id','').strip()
        record_name = self.get_argument('record_name','').strip()
        record_description = self.get_argument('record_description','').strip()
        digest_enzyme = self.get_argument('digest_enzyme','').strip()
        
        ##Randomly assign name for temporary directory
        tmp_str = base64.b32encode(os.urandom(20)).decode('utf-8')
        tmp_dir_name = "./output_folders/"+tmp_str
        if not os.path.isdir("./output_folders/"):
            os.mkdir("./output_folders/")
        os.mkdir(tmp_dir_name)
        current_dir = os.getcwd()
    
        ##Parse text box
        if len(text_box) == 0:
            self.finish("No MoClo part string provided")
            return
        else:
            text_box = text_box.strip()
            text_box = text_box.replace("\r","\n")
            text_split = text_box.split("\n")
            text_split[:] = [item for item in text_split if item != ''] ##remove empty strings
    
        self.write("Loaded "+str(len(text_split))+" record(s) successfully<br>")        
        plasmidNames = auto_golden_gate_functions.reducedPlasmidNamesToCanonical(text_split)
        self.write("Parsed input:<br>")
        self.write(str(plasmidNames)+"<br>")

        os.chdir(tmp_dir_name)
        files = glob.glob("../../pydna_development/*.gb*")
        self.write(str(len(files))+" files found.<br>")
 
        plasmidSeqs,featuresToApply =  auto_golden_gate_functions.extractPlasmidSeqsAndFeatures(plasmidNames,files)
        self.write("Loaded "+str(len(plasmidSeqs))+" plasmid files with "+str(len(featuresToApply))+" annotations.<br>")
        if len(plasmidNames) != len(plasmidSeqs):
            self.write("Error: Number of plasmids parsed is less than requested.<br>")
            self.write("This means the program couldn't find some file(s).<br>")
            self.finish("exiting...<br>")
            return
        
        self.write("Performing in-silico "+digest_enzyme+" digestion...<br>")
        fragments =  auto_golden_gate_functions.digestBsaI(plasmidSeqs,digest_enzyme) ##Output in a dictonary with some additional information
        initial_frags = [item for sublist in fragments.values() for item in sublist] ##Get rid of the dictionary structure
        self.write("Digestion resulted in "+str(len(initial_frags))+ " fragments.<br>")
        initial_frags_filtered =  auto_golden_gate_functions.subtract_frags_by_seq(initial_frags) ##Subtract away GFP and chloramphenicol resistance

        pastFrags = []
        keepLooping = True
        i=0
        while keepLooping == True:
            i+=1
            self.write("Ligation loop:"+str(i)+"<br>")
            combinedFrags = auto_golden_gate_functions.ligate_fragments(pastFrags+initial_frags_filtered)
            self.write("Number of resulting fragments:"+str(len(combinedFrags))+"<br>")

            for f in combinedFrags:
                if f.circular:
                    keepLooping = False
                    self.write("Found circular plasmid. Exiting ligation loop.<br>")
                    break
            pastFrags = combinedFrags
            if i > 19:
                os.chdir(current_dir)
                self.write("Reached limit of 20 ligation cycles. Seems ligation is not converging to circular plasmids.<br>")
                self.finish("Exiting...<br>")
                return 
        self.write("Number of ligation products:"+str(len(combinedFrags))+"<br>")
        circularPlasmids = auto_golden_gate_functions.filterToCircularDSeqs(combinedFrags)
        self.write("Number of ligation products after fitering to circular only:"+str(len(circularPlasmids))+"<br>")
        dereplicatedPlasmids = auto_golden_gate_functions.dereplicate_circular_sequences(circularPlasmids)
        self.write("Number of ligation products after dereplicating:"+str(len(dereplicatedPlasmids))+"<br>")
        assert len(dereplicatedPlasmids) == 1

        self.write("Reapplying feature annotations to sequence and exporting...<br>")
        from Bio.Alphabet.IUPAC     import IUPACAmbiguousDNA
        import datetime
        theSeq = Bio.Seq.Seq(str(dereplicatedPlasmids[0]),IUPACAmbiguousDNA()).upper()
        theRecord = Bio.SeqRecord.SeqRecord(theSeq,id="testID", name="testName",description="testDescription")
        theRecord.annotations["topology"] = "circular"
        theRecord.annotations["molecule_type"] = "ds-DNA"
        today = datetime.datetime.today()
        today_string = today.strftime("%d-%b-%Y").upper() ##%d-%b-%Y -> uppercase
        theRecord.annotations["date"] = today_string

        theRecord = auto_golden_gate_functions.reapply_features_to_record(theRecord,featuresToApply)

        Bio.SeqIO.write(theRecord, record_name+".gbk", "gb")

        os.chdir(current_dir)
        
        self.write("<br>Produced plasmid file(s):")
        zipname = tmp_dir_name+"/plasmids.zip"
        zip_handle = None
        for root, subFolders, files in os.walk(tmp_dir_name):
            for filename in files:
                if filename.lower().endswith(".gbk"):
                    self.write("<br>"+filename)
                    if zip_handle == None: 
                        zip_handle = zipfile.ZipFile(zipname,"w") ##Open the zip file.
                    zip_handle.write(tmp_dir_name+"/"+filename,filename) ##Append files to the zip file.
        
        if not zip_handle == None:
            zip_handle.close()
        
        if os.path.isfile(zipname):
        ##Check if the file was actually produced
            self.write("<br><a href="+"output_folders/"+tmp_str+"/plasmids.zip>Plasmid .zip download</a>")
        else:
            self.write("<br>Couldn't produce plasmid files (likely due to missing reference .gb internally)")
    
    
        self.write("</div>")

        ##Making the return text.  Tries to return only unique records (non-matching SeqID & sequence)
        ##Since it holds the sequences in memory to compare if duplicates exist, this could deal badly with large memory fasta files.
        self.finish("</font></html>")
        #self.finish(cname + " is uploaded!! Check %s folder" %__UPLOADS__)
        ##return 0
 
 
application = tornado.web.Application([
        (r"/", Userform),
        (r"/upload", Upload),
        (r"/output_folders/(.*)",tornado.web.StaticFileHandler,{'path':"./output_folders/"}) ##http://stackoverflow.com/questions/10165665/using-tornado-how-do-i-serve-static-files-and-serve-a-favicon-ico-from-a-differ
        ], debug=True)
 
 
if __name__ == "__main__":
    
    pid = str(os.getpid())
    f = open('fasta_lookup_server.pid', 'w')
    f.write(pid)
    f.close()
    print("Running server on TCP_PORT:",TCP_PORT)
    application.listen(TCP_PORT)
    tornado.ioloop.IOLoop.instance().start()
