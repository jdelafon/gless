import Tkinter as tk
import os,sys
import argparse,re
#fuckinfile = "/Users/julien/Library/Saved Application State/org.python.python.savedState"
#os.chmod(fuckinfile,0o777)
#os.remove(fuckinfile) # bug on osx

class Parse:
    """A replacement for btrack when bbcflib is not found, able to parse
       bed and bedGraph formats only."""
    def __init__(self,filename):
        self.filehandle = os.path.abspath(filename)
        self.format = os.path.splitext(filename)[1]
    def read(self,chrom=None,fields=None):
        with open(self.filehandle) as f:
            for line in f:
                line = line.strip().split()
                try:
                    chr = line[0]
                    start,end = (int(line[1]),int(line[2]))
                except IndexError:
                    raise ValueError(("Library 'bbcflib' not found. "
                                      "Only 'bed' and 'bedGraph' formats available "
                                      "(got '%s')." % os.path.basename(self.filehandle)))
                if chrom and chr != chrom:
                    continue
                yield (start,end)+tuple(line[3:])
try:
    from bbcflib.btrack import track
    assert track
except ImportError:
    track = Parse


def read10(trackList, nfeat=None, nbp=None, selection=None):
    """Yield a list of lists [[(1,2),(3,4)], [(1,3),(5,6)]] with either the *nfeat*
       next items, or all next items within an *nbp* window."""
    chr = 'chr1'
    fields = ['start','end']
    streams = [track(t).read(chr,fields=fields) for t in trackList]
    available_streams = range(len(streams))
    if nbp: nfeat = None
    elif nfeat: nbp = None
    if nfeat:
        current = [(s.next(),k) for k,s in enumerate(streams)]
        sortkey = lambda x:x[0][0]
        # Repeat each time the function is called
        while available_streams:
            toyield = [[] for _ in streams]
            len_toyield = 0
            # Yield *nfeat* items
            while len_toyield < nfeat and current:
                current.sort(key=sortkey)
                min_idx = current[0][-1]
                toyield[min_idx].append(current.pop(0)[0])
                len_toyield += 1
                try: current.append([streams[min_idx].next(),min_idx])
                except StopIteration:
                    try: available_streams.pop(min_idx)
                    except IndexError: continue
            # Add feats that go partially beyond
            if any(toyield):
                maxpos = max(x[-1][1] for x in toyield if x)
                for n,x in enumerate(current):
                    if x[0][0] < maxpos:
                        toyield[x[-1]].append((x[0][0],min(x[0][1],maxpos)))
                        current[n][0] = (min(x[0][1],maxpos),x[0][1])
                yield toyield
            else: break
    elif nbp:
        current = [s.next() for k,s in enumerate(streams)]
        sortkey = lambda x:x[0][0]
        # Repeat each time the function is called
        k = 0 # number of times called
        while available_streams:
            k += 1
            toyield = [[] for _ in streams]
            toremove = []
            # Load items within *nbp* in *toyield*
            for n in available_streams:
                while current[n]:
                    x = current[n]
                    if x[1] <= k*nbp:
                        toyield[n].append((x[0]-(k-1)*nbp,x[1]-(k-1)*nbp))
                        try: current[n] = streams[n].next()
                        except StopIteration:
                            current[n] = None
                            toremove.append(n)
                    elif x[0] < k*nbp:
                        toyield[n].append((x[0]-(k-1)*nbp,nbp))
                        current[n] = (k*nbp,x[1])
                    else: break
            for n in toremove: available_streams.remove(n)
            if any(toyield):
                yield toyield
            else:
                yield [[(0,0)] for _ in streams]

def draw(tracks_content, geometry=None, WIN_WIDTH=700, nfeat=None, nbp=None):
    """Create a new window and draw from the *tracks_content* coordinates
       (of the form [[(1,2),(3,4)], [(3,5),(4,6)],...] )."""
    def _bp2px(y,wwidth,reg_bp):
        try: return y*wwidth/reg_bp
        except ZeroDivisionError: return 0
    def exit_on_Esc(event):
        if event.keysym == 'Escape':
            root.destroy()
            sys.exit(0)
    def shift_on_arrow(event):
        if event.keysym == 'Left':
            pass # backwards??
        if event.keysym == 'Right':
            pass
    root = tk.Tk()
    root.bind("<Key>", exit_on_Esc)
    if geometry: # keep previous position of the window
        geometry = '+'.join(['']+geometry.split('+')[1:]) # "+70+27"
        root.geometry(geometry)
    htrack = 30
    feat_pad = 10
    feat_len = htrack - 2*feat_pad
    reg_bp = max(x[1] for t in tracks_content for x in t) # whole region size in bp
    reg_bp = max(reg_bp,nbp)
    for n,t in enumerate(tracks_content):
        l = tk.Label(root,text="track%d"%(n+1))
        lwidth = l.winfo_reqwidth()
        l.grid(row=n,column=0,padx=5,pady=5)
        c = tk.Canvas(root, width=WIN_WIDTH-lwidth, height=htrack, bg='white')
        c.grid(row=n,column=1)
        for k,feat in enumerate(t):
            x1,x2 = (feat[0],feat[1])
            y1,y2 = (0+feat_pad,feat_len+feat_pad)
            x1 = _bp2px(x1,WIN_WIDTH-lwidth,reg_bp)
            x2 = _bp2px(x2,WIN_WIDTH-lwidth,reg_bp)
            c.create_rectangle(x1+1,y1,x2+1,y2, fill="blue")
    root.wm_attributes("-topmost", 1) # makes the window appear on top
    return root


def gless(trackList, nfeat=6, nbp=None, selection=None):
    """Main controller function after option parsing."""
    tracks = read10(trackList,nfeat=nfeat,nbp=nbp,selection=selection)
    try: tracks_content = tracks.next()
    except StopIteration:
        print "Nothing to show"
        sys.exit(0)
    drawn = 0
    geometry = None
    while True:
        if not drawn:
            root = draw(tracks_content,geometry,nfeat=nfeat,nbp=nbp)
            geometry = root.geometry()
        drawn += 1
        key = raw_input()
        if key == '': sys.exit(0) # "Enter" pressed
        elif key == ' ':
            try:
                tracks_content = tracks.next()
                root.destroy()
                drawn = 0
            except StopIteration: print "End of file."

def main():
    parser = argparse.ArgumentParser(description="Graphical 'less' for track files.")
    parser.add_argument('-n','--nfeat', default=None, type=int,
                       help="Number of features to display.")
    parser.add_argument('-b','--nbp', default=None, type=int,
                       help="Number of base pairs to display.")
    parser.add_argument('-s','--sel', default=None,
                       help="Region to display, formatted as in 'chr1:12-34'.")
    parser.add_argument('file', nargs='+', default=None,
                       help='A set of track files, separated by spaces')
    args = parser.parse_args()
    trackList = args.file
    trackList = ['testing_files/test1.bed','testing_files/test2.bed']
    assert (bool(args.nfeat) != bool(args.nbp)) != bool(args.sel) # one at a time (XOR)
    if args.sel: assert re.search('chr[0-9XY]+:[0-9XY]+-[0-9XY]+',args.sel)
    gless(trackList,nfeat=args.nfeat,nbp=args.nbp,selection=args.sel)

if __name__ == '__main__':
    sys.exit(main())


#root.focus_set()
#print c.winfo_reqheight(), c.winfo_reqwidth()
#print c.winfo_width(), c.winfo_height()
#root.mainloop()
#print 'Toyield',str(toyield)
#print 'Toyield', [[(x[0]+(k-1)*nbp,x[1]+(k-1)*nbp) for x in t] for t in toyield]
#print "Window", str((k-1)*nbp)+'-'+str(k*nbp)
