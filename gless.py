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
        self.format = os.path.splitext(filename)[1][1:]
    def read(self,fields=None,selection=None):
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
                if selection and selection['chr'] != chr:
                    continue
                yield (start,end)+tuple(line[3:])
try:
    from bbcflib.btrack import track
    assert track
except ImportError:
    track = Parse

def format(filename):
    """Return whether it is a track with 'intervals' or a 'density'."""
    with track(filename) as t:
        if t.format.lower() in ['bed']:
            return 'intervals'
        elif t.format.lower() in ['bedgraph','wig','bigWig','sga']:
            return 'density'

def parse_selection(sel):
    if not sel: return
    elif re.search('^chr[0-9XY]*:[0-9XY]+-[0-9XY]+',sel):
        chr,coord = sel.split(':')
        st,en = coord.split('-')
        return {'chr':chr,'start':(int(st),int(en)),'end':(int(st),int(en))}
    elif re.search('^chr[0-9XY]*$',sel):
        return {'chr':sel}
    else: print "Bad region formatting, got -s %s,." % sel; sys.exit(1)

def read10(trackList,nfeat,nbp,sel):
    """Yield a list of lists [[(1,2,n),(3,4,n)], [(1,3,n),(5,6),n]] with either the *nfeat*
       next items, or all next items within an *nbp* window. `n` is a name or a score."""
    #sel = 'chr1'
    tracks = [track(t) for t in trackList]
    sel = parse_selection(sel)
    streams = [t.read(fields=t.fields[1:4],selection=sel) for t in tracks]
    available_streams = range(len(streams))
    if nbp: nfeat = None
    elif nfeat: nbp = None
    if nfeat:
        current = [[s.next(),k] for k,s in enumerate(streams)]
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
                        toyield[x[-1]].append((x[0][0],min(x[0][1],maxpos),x[0][2]))
                        current[n][0] = [min(x[0][1],maxpos),x[0][1],x[0][2]]
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
                        toyield[n].append((x[0]-(k-1)*nbp,x[1]-(k-1)*nbp,x[2]))
                        try: current[n] = streams[n].next()
                        except StopIteration:
                            current[n] = None
                            toremove.append(n)
                    elif x[0] < k*nbp:
                        toyield[n].append((x[0]-(k-1)*nbp,nbp,x[2]))
                        current[n] = (k*nbp,x[1],x[2])
                    else: break
            for n in toremove: available_streams.remove(n)
            if any(toyield):
                yield toyield
            else:
                yield [[(0,0,'')] for _ in streams]

def draw(names,tracks_content,geometry,nfeat,nbp,ntimes,maxpos):
    """Create a new window and draw from the *tracks_content* coordinates
       (of the form [[(1,2,n),(3,4,n)], [(3,5,n),(4,6,n)],...],
       where `n` is either a name or a score)."""
    WIN_WIDTH = 800
    htrack = 30
    rmargin = 100
    feat_pad = 10
    bg = "grey"
    canvas_bg = "white"
    line_col = "black"
    feat_col = "blue"
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
    root.config(bg=bg)
    if geometry: # keep previous position of the window
        geometry = '+'.join(['']+geometry.split('+')[1:]) # "+70+27"
        root.geometry(geometry)
    feat_thk = htrack - 2*feat_pad
    reg_bp = max(x[1] for t in tracks_content for x in t) # whole region size in bp
    reg_bp = max(reg_bp,nbp)
    canvas = [tk.Canvas(root,height=htrack,bd=0,bg=canvas_bg,highlightthickness=0)
              for _ in names]
    wlabel = 0
    # Draw labels
    for n,name in enumerate(names):
        l = tk.Label(root,text=names[n],bd=0,highlightthickness=0,bg=bg,padx=5)
        wlabel = max(l.winfo_reqwidth(),wlabel)
        l.grid(row=n,column=0)#,ipadx=4,ipady=4)
    wcanvas = WIN_WIDTH-wlabel-rmargin
    # Draw the tracks
    for n,t in enumerate(tracks_content):
        type = format(names[n])
        c = canvas[n]
        c.config(width=wcanvas)
        c.grid(row=n,column=1)
        if type == 'intervals':
            c.create_line(0,htrack/2.,WIN_WIDTH,htrack/2.,fill=line_col) # track axis
            y1,y2 = (0+feat_pad,feat_thk+feat_pad)
            for k,feat in enumerate(t):
                x1,x2,g = (feat[0],feat[1],feat[2])
                x1 = _bp2px(x1,wcanvas,reg_bp)
                x2 = _bp2px(x2,wcanvas,reg_bp)
                c.create_rectangle(x1,y1,x2,y2, fill=feat_col)
        elif type == 'density':
            if t: top_bp = max(float(x[2]) for x in t) # highest score
            c.create_line(0,htrack-1,WIN_WIDTH,htrack-1,fill=line_col) # track axis
            for k,feat in enumerate(t):
                x1,x2,s = (feat[0],feat[1],feat[2])
                x1 = _bp2px(x1,wcanvas,reg_bp)
                x2 = _bp2px(x2,wcanvas,reg_bp)
                s = _bp2px(s,htrack,top_bp)
                c.create_rectangle(x1,htrack-1,x2,htrack-s+5,fill=feat_col)
    # Add a blank frame on the right as a margin
    for n in range(len(names)):
        w = tk.Frame(width=rmargin,height=htrack,bg=bg)
        w.grid(row=n,column=2)
    # Axis
    minpos = maxpos
    if nbp: maxpos = ntimes*nbp
    elif nfeat: maxpos = reg_bp
    min_label = tk.Label(text=str(minpos),bd=0,bg=bg,anchor='e')
    min_label.grid(row=n+1,column=0,sticky='e',padx=5)
    c = tk.Canvas(root,width=wcanvas,height=2*htrack,bd=0,bg=canvas_bg,highlightthickness=0)
    c.grid(row=n+1,column=1)
    pad = c.winfo_reqheight()/2.
    c.create_line(0,pad,WIN_WIDTH,pad,fill=line_col)  # axis
    # Ticks & labels
    for n,t in enumerate(tracks_content):
        for k,feat in enumerate(t):
            f1,f2 = (feat[0],feat[1])
            x1 = _bp2px(f1,wcanvas,reg_bp)
            x2 = _bp2px(f2,wcanvas,reg_bp)
            c.create_line(x1,pad,x1,pad-5,fill=line_col)
            c.create_line(x2,pad,x2,pad+5,fill=line_col)
            if nbp:
                f1 = feat[0]+(ntimes-1)*nbp
                f2 = feat[1]+(ntimes-1)*nbp
            if f1!=minpos and f1!=maxpos:
                c.create_text(x1,pad-5,text=str(f1),anchor='s')
            if f2!=minpos and f2!=maxpos:
                c.create_text(x2,pad+5,text=str(f2),anchor='n')
    max_label = tk.Label(text=str(maxpos),bd=0,bg=bg,anchor='w')
    max_label.grid(row=n+1,column=2,sticky='w',padx=5)
    root.wm_attributes("-topmost", 1) # makes the window appear on top
    return root,maxpos

def gless(trackList, nfeat=None, nbp=None, sel=None):
    """Main controller function after option parsing."""
    names = [os.path.basename(t) for t in trackList]
    tracks = read10(trackList,nfeat,nbp,sel)
    try: tracks_content = tracks.next()
    except StopIteration:
        print "Nothing to show"
        sys.exit(0)
    needtodraw = True
    ntimes = 1 # Number of times spacebar is hit
    maxpos = 0 # Last biggest coordinate
    geometry = None
    while True:
        if needtodraw:
            root,maxpos = draw(names,tracks_content,geometry,nfeat,nbp,ntimes,maxpos)
            geometry = root.geometry()
            needtodraw = False
            ntimes += 1
        key = raw_input()
        if key == '': sys.exit(0) # "Enter" pressed
        elif key == ' ':
            try:
                tracks_content = tracks.next()
                root.destroy()
                needtodraw = True
            except StopIteration: print "End of file."

def main():
    parser = argparse.ArgumentParser(description="Graphical 'less' for track files.")
    parser.add_argument('-n','--nfeat', default=10, type=int,
                       help="Number of features to display.")
    parser.add_argument('-b','--nbp', default=None, type=int,
                       help="Number of base pairs to display.")
    parser.add_argument('-s','--sel', default=None,
                       help="Region to display, formatted as in 'chr1:12-34'.")
    parser.add_argument('file', nargs='+', default=None,
                       help='A set of track files, separated by spaces')
    args = parser.parse_args()
    if args.nfeat and args.nbp:
        print "Only one of -n/-b is allowed at a time."; sys.exit(1)
    elif args.nbp: args.nfeat = None
    elif args.nfeat and args.nfeat > 10000:
        print "Up to 10000 features permitted, got -n %s." % args.nfeat; sys.exit(1)
    gless(trackList=args.file,nfeat=args.nfeat,nbp=args.nbp,sel=args.sel)

if __name__ == '__main__':
    sys.exit(main())

# python gless.py -n 12 testing_files/test1.bed testing_files/test2.bed testing_files/yeast_genes.bed


#trackList = ['testing_files/test1.bed','testing_files/test2.bed']
#root.focus_set()
#print c.winfo_reqheight(), c.winfo_reqwidth()
#print c.winfo_width(), c.winfo_height()
#root.mainloop()
#print 'Toyield',str(toyield)
#print 'Toyield', [[(x[0]+(k-1)*nbp,x[1]+(k-1)*nbp) for x in t] for t in toyield]
#bg = "#%02x%02x%02x" % (255, 255, 224) #beige background
#print "Window", str((k-1)*nbp)+'-'+str(k*nbp)
        #c.create_line(0,0,0,htrack,fill="grey") # separator label|canvas
    #c.create_line(0,0,0,2*pad,fill=line_col)         # separator
