    
class Fasta:
    """"This is a class made to hold multiple contigs."""
    def __init__(self, contigs = [], species = None):
        self.contigs = contigs # this is a list of Contig objects
        self.species = species # this is metadata in case we want it
        self.length = len(contigs)
    
    def __str__(self):
      return f"Fasta species: {self.species}\nFasta length: {self.length}."
        
    def dump(self):
        # When used with print(), it dumps an entire formatted name and seq fasta combo to stdout.
        # This is useful for outputting the entire fasta without having to use a list (like with Contig class)
        for contig in self.contigs:
            print(contig.dump())

class Contig:
    """This is a class for holding contig data/names from fasta files. Yes, this is a worse reinvention of SeqIO, but it's also more portable."""
    def __init__(self, name = None, seq = None, length = 0):
        self.name = name
        self.seq = seq
        self.length = length

    def __str__(self):
      return f"Contig name: {self.name}\nContig length: {self.length}."
  
    def to_dict(self):
        return {
            'name': self.name,
            'seq': self.seq,
            'length': self.length
            }
    
    def reverseComplement(self):
        # reverse complements your string
        complement_dict = {"A": "T", "a": "t", # four standard nucleotides
                           "T": "A", "t": "a", 
                           "C": "G", "c": "g",
                           "G": "C", "g": "c",
                           # I don't think I'll have many ambiguous bases, but...
                           "R": "Y", "r": "y",
                           "Y": "R", "y": "r",
                           "K": "M", "k": "m",
                           "M": "K", "m": "k",
                           "B": "V", "b": "v",
                           "V": "B", "v": "b",
                           "D": "H", "d": "h",
                           "H": "D", "h": "d"
                           }
        complement_trans_tbl = self.seq.maketrans(complement_dict)
        self.seq = self.seq.translate(complement_trans_tbl)[::-1] # inverts and complements string
    
    def dump(self):
        # When used with print(), it dumps an entire formatted name and seq fasta combo to stdout.
        # This is useful for outputting the fasta.
        return f">{self.name}\n{self.seq}"
    

def parse_fasta(fasta): 
    ctg_lst= [] # initialize a list 
    seq = name = ""
    
    with open(fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if seq == "":
                    name = line.strip().lstrip('>')
                else:
                    #initialize new contig and write to list
                    Ctg = Contig(name = name, seq = seq, length = len(seq))
                    ctg_lst.append(Ctg)
                    name = seq = "" #reset the loop for the next round
                    name = line.strip().lstrip('>')
            else:
                seq = seq + line.strip()
        Ctg = Contig(name = name, seq = seq, length = len(seq))
        ctg_lst.append(Ctg)
    return ctg_lst
    
def sortSeq_lst(seq):
    # Small function that provides a key to sort list objects
    # Works for Alignment and Contig objects because they both have lengths
    return seq.length

