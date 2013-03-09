import Tkinter as tk
import os,sys

class FakeTrack:
    """A replacement for btrack when bbcflib is not found, able to parse
       bed and bedGraph formats only."""
    def __init__(self,filename):
        self.filehandle = filename
        self.format = os.path.splitext(filename)[1]
    def read(self,chrom=None,fields=None):
        with open(self.filehandle) as f:
            for line in f:
                line = line.strip().split()
                chr = line[0]
                start,end = (int(line[1]),int(line[2]))
                if chrom and chr != chrom:
                    continue
                yield (start,end)+tuple(line[3:])

try:
    from bbcflib.btrack import track
    assert track
except ImportError:
    print "Library 'bbcflib' not found. Only .bed and .bedGraph formats available."
    track = FakeTrack


def read_tracks(trackList, nfeat=None, nbp=None):
    """Yield a list of lists [[(1,2),(3,4)], [(1,3),(5,6)]] with either the *nfeat*
       next items, or all next items within an *nbp* window."""
    chr = 'chr1'
    fields = ['start','end']
    streams = [track(t).read(chr,fields=fields) for t in trackList]
    available_streams = range(len(streams))
    current = [(s.next(),k) for k,s in enumerate(streams)]
    if nfeat:
        # Repeat each time the function is called
        sortkey=lambda x:x[0][0]
        while available_streams:
            toyield = [[] for _ in streams]
            len_toyield = 0
            # Yield *nfeat* items
            while len_toyield < nfeat and current:
                current.sort(key=sortkey)
                #print 'current:',current
                min_idx = current[0][-1]
                toyield[min_idx].append(current.pop(0)[0])
                len_toyield += 1
                try:
                    current.append([streams[min_idx].next(),min_idx])
                except StopIteration:
                    try:
                        available_streams.pop(min_idx)
                    except IndexError:
                        continue
                #print 'toyield',toyield
            # Add feats that go partially beyond
            if any(toyield):
                maxpos = max(x[-1][1] for x in toyield if x)
                for n,x in enumerate(current):
                    if x[0][0] < maxpos:
                        toyield[x[-1]].append((x[0][0],min(x[0][1],maxpos)))
                        current[n][0] = (min(x[0][1],maxpos),x[0][1])
                #print "Yielded",toyield
                yield toyield
            else: break
    elif nbp:
        pass

def pos2px(y,wwidth,reg_bp):
    return y*wwidth/reg_bp

def draw(tracks_content,geometry=None):
    def exit_on_keydown(event):
        if event.char == chr(27):
            sys.exit(0)
    root = tk.Tk()
    if geometry: # keep previous position of the window
        geometry = '+'.join(['']+geometry.split('+')[1:]) # "+70+27"
        root.geometry(geometry)
    root.bind("<Key>", exit_on_keydown)
    wwindow = 700
    htrack = 30
    feat_pad = 10
    feat_len = htrack - 2*feat_pad
    reg_bp = max(x[1] for t in tracks_content for x in t) # whole region size in bp
    for n,t in enumerate(tracks_content):
        l = tk.Label(root,text="track%d"%(n+1))
        lwidth = l.winfo_reqwidth()
        l.grid(row=n,column=0,padx=5,pady=5)
        c = tk.Canvas(root, width=wwindow-lwidth, height=htrack, bg='white')
        c.grid(row=n,column=1)
        for k,feat in enumerate(t):
            x1,x2 = (feat[0]+1,feat[1]+1)
            y1,y2 = (0+feat_pad,feat_len+feat_pad)
            x1 = pos2px(x1,wwindow-lwidth,reg_bp)
            x2 = pos2px(x2,wwindow-lwidth,reg_bp)
            c.create_rectangle(x1,y1,x2,y2, fill="blue")
    root.wm_attributes("-topmost", 1) # makes the window appear on top
    #root.mainloop()
    return root

def gless(trackList, nfeat=6, nbp=None):
    tracks = read_tracks(trackList,nfeat=nfeat,nbp=nbp)
    try:
        tracks_content = tracks.next()
    except StopIteration:
        print "Nothing to show"
        sys.exit(0)
    drawn = 0
    geometry = None
    while True:
        if not drawn:
            root = draw(tracks_content,geometry)
            geometry = root.geometry()
        drawn += 1
        key = raw_input()
        if key == chr(27):
            sys.exit(0)
        elif key == ' ':
            try:
                tracks_content = tracks.next()
                root.destroy()
                drawn = 0
            except StopIteration:
                print "Nothing more to display."
                sys.exit(0)

trackList = ['testing_files/test1.bed','testing_files/test2.bed']
gless(trackList)


        #root.focus_set()
        #print c.winfo_reqheight(), c.winfo_reqwidth()
        #print c.winfo_width(), c.winfo_height()

#tracks = [[('chr1',1,12,6.0),('chr1',19,34,2.0),('chr1',50,52,9.0)],
#          [('chr1',11,17),('chr1',23,75)],
#         ]
#graph(tracks)
