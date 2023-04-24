class Protein:
    def __init__(self, pdbid, chain, knot, knot_type):
        self.pdbid = pdbid
        self.chain = chain
        self.knot = knot
        self.knot_type = knot_type

    def add_data(self, data):
        self.data = data

    def add_barcode(self, barcode):
        self.barcode = barcode

    def add_pl(self, pl):
        self.pl = pl