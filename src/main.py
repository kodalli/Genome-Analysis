"""
    Genome Analysis Project
"""

from string import digits
from Bio.Seq import Seq


# Hepatitis delta virus isolate GZ37, complete genome

hepdRaw = """1 catgagccac catccgaacg aagattgcgc gaggggcggg atcagcgccc gagaggggta
       61 agtggtaaag agcattggaa cgtcggagaa actactccca agaaggaaaa aagagaaagc
      121 aagaatcgga cgagttcccc aagacgctgg gaacgtctcg gaaggggaaa gaaggaaggt
      181 gggaaagaaa ggggcgggcc tcccgatccg aggggcccaa ccaccaagtt tggagagcac
      241 tccggcccca agggttgaga gtacccagaa ggaggaatcc tctcggagaa aagcagataa
      301 atcacctcca gaggacccct tcagcgaaca aaggagctct gaagcgcgag gagtaagagc
      361 atagcgatag ggggagatgc taggagttag gggagaccga agcgaggagg aaagcaaaga
      421 aagcaacggg gctagccggt aggtgttccg cctcccgaga ggggacgagt gaggcttatc
      481 ccggggaact cggcgaatcg ttcccacata gcagacccca ggaccccctt ccaaatggtc
      541 cgaggggggt ggctaggaac acaggggacc ggtggagcca tgggatgctc ctcccgatgt
      601 ccgaatccat ccctcccccc gagggtcgcc caggaatggc gggaccccac tcaactgggg
      661 tccgcgttcc atcctttctt acctgatggc cggcatggtc ccagcctcct cgctggcgcc
      721 ggctgggcaa cattccgagg ggaccgtccc tcggtaatgg cgaatgggac ccagaaatct
      781 ctctagattc ccagagagaa tcgagagaaa actggctctc ccttagccat ccgagtggac
      841 gttcgtcctc cttcggatgc ccaggtcgga ccgcgaggag gtggagatgc catgccgacc
      901 cgaagaggaa agaaggacgc gagacacgaa cccgtgagtg gaaacccgct ttattcactg
      961 gggtcgacaa ctctggggag aaaagggagg atcggctggg aagagtatat cctatgggaa
     1021 tccccggtct ccccttatgt ccagcccctc cccggtcctg gtgaaggggg actccggaat
     1081 tccttgcatg ccgggaacga agccgccccc gggcgctccc ctcgatccac cttcgagggg
     1141 gttcacacct ccaaccggcg ggccggctac tcttctttcc cttctctcgt cttcctcggt
     1201 caacctctta agttcctctt cttcctcctt gctgagcttc ttccctccgg cactcagctg
     1261 cttccttttg ttctcgaggg ccttccttcg tcggtgatcc tgcctctcct tgtcggagaa
     1321 ccctcctctg agaggcctct tcctaggtcc ggagtctacc tccatctggt ctgttcgggc
     1381 cctcttcgcc gggggagccc cctctccatc cttatccttc tttccgagaa ttcctttgat
     1441 gtttcccagc cagggatttt cgtcctcaag tttcttgatt ttcttcttaa ccttccggag
     1501 gtctctctcg agttcctcta acttctttct tccactcacc cactgctcga gaacctcttc
     1561 cctgcccccg cgatttttct tcgattcgga gcggctcatc tcgacaagag gcggggctcc
     1621 ctagttctct tactcttttc tgtaaagagg agactgtggg agtccccgcc caagttcgag """

hepdParsed  = hepdRaw.translate(str.maketrans('', '', digits)).replace(' ', '').upper()
hepd_seq = Seq(hepdParsed)

print(hepd_seq.count("A"))
print(hepd_seq.complement())