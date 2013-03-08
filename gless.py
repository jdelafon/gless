import Tkinter as tk
import os,sys
from operator import itemgetter

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
    import iuqwyer
except ImportError:
    print "Library 'bbcflib' not found. Only .bed and .bedGraph formats available."
    track = FakeTrack


def read_tracks(trackList, nfeat=10, nbp=None):
    """Yield a list of lists [[(1,2),(3,4)], [(1,3),(5,6)]] with either the *nfeat*
       next items, or all next items within an *nbp* window."""
    if nbp: nfeat = sys.maxint
    elif nfeat: nbp = sys.maxint
    sortkey=lambda x:x[0][0]
    chr = 'chr1'
    fields = ['start','end']
    streams = [track(t).read(chr,fields=fields) for t in trackList]
    available_streams = range(len(streams))
    current = [(s.next(),k) for k,s in enumerate(streams)]
    current.sort(key=sortkey)
    if nfeat:
        # Repeat each time the function is called
        while available_streams:
            toyield = [[] for _ in streams]
            # Yield *nfeat* items
            while len(toyield) < nfeat and current:
                #print 'current:',current
                min_idx = current[0][-1]
                toyield[min_idx].append(current.pop(0)[0])
                try:
                    current.append((streams[min_idx].next(),min_idx))
                except StopIteration:
                    try:
                        available_streams.pop(min_idx)
                    except IndexError:
                        continue
                current.sort(key=sortkey)
                #print "toyield:",toyield
            yield toyield

def pos2px(y,wwidth,reg_bp):
    return y*wwidth/reg_bp

def graph(tracks, wwidth=375, wheight=325):
    root = tk.Tk()
    leftm = 50
    rightm = 5
    topm = 50
    botm = 50
    track_width = 30
    feat_pad = 10
    feat_width = track_width - 2*feat_pad
    reg_bp = max(x[1] for t in tracks for x in t) # whole region size in bp
    for n,t in enumerate(tracks):
        l = tk.Label(root,text="track%d"%(n+1))
        l.grid(row=n,column=0,padx=5,pady=5)
        c = tk.Canvas(root, height=track_width, bg='white')
        print c.winfo_reqheight(), c.winfo_reqwidth()
        print c.winfo_width(), c.winfo_height()
        c.grid(row=n,column=1)
        for k,feat in enumerate(t):
            x1,x2 = (feat[0]+1,feat[1]+1)
            y1,y2 = (0+feat_pad,feat_width+feat_pad)
            x1 = pos2px(x1,wwidth,reg_bp)
            x2 = pos2px(x2,wwidth,reg_bp)
            c.create_rectangle(x1,y1,x2,y2, fill="blue")
    root.mainloop()


trackList = ['testing_files/test1.bed','testing_files/test2.bed']
res = read_tracks(trackList,nfeat=10)
r = res.next()
print r
graph(r)



#tracks = [[('chr1',1,12,6.0),('chr1',19,34,2.0),('chr1',50,52,9.0)],
#          [('chr1',11,17),('chr1',23,75)],
#         ]
#graph(tracks)
