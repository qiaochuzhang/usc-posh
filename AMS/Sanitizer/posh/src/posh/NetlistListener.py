from netlist_parser.SpectreParserListener import SpectreParserListener


class NetlistListener(SpectreParserListener):
    def enterNetlist(self, ctx):
        print "Enter netlist"

    def exitNetlist(self, ctx):
        print "Exit netlist"

    def enterEveryRule(self, ctx):
        #print "All", type(ctx)
        pass