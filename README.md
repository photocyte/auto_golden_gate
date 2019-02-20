# auto_golden_gate
A Python based server for making Golden Gate cloning a little faster

This is a standalone python based server that intends to make the production of [Gibson assembly](https://en.wikipedia.org/wiki/Gibson_assembly) primers & in-silico plasmids from that cloning more straightforward. In-silico plasmids are produced in GenBank (.gb) format with an eye toward proper integration with [SnapGene](http://www.snapgene.com/).  Appropriate insert annotations (assuming the sequences you input are [CDSs](https://en.wikipedia.org/wiki/Coding_region)) are produced in the resulting .gb file(s). Primers are output either in CSV format suitable for copy pasting, or 96-well plate tables in a Excel .xlsx, both formats suitable for direct oligonucleotide synthesis orders from [IDT](https://www.idtdna.com/pages).

### Dependencies:
 * Python3.6+
 * BioPython
 * SciPy (for the Tornado stand alone webserver)
 * [pydna](https://github.com/BjornFJohansson/pydna)

