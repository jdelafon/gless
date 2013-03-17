import Tkinter as tk
import os,sys
import argparse,re

###############################################################################

class Parser(object):
    """A replacement for btrack when bbcflib is not found, able to parse
       bed and bedGraph formats only."""
    def __init__(self,filename):
        self.path = os.path.abspath(filename)
        self.format = os.path.splitext(filename)[1][1:]
        self.fields = ['chr','start','end','other']
    def __enter__(self):
        return self
    def __exit__(self,errtype,value,traceback):
        pass
    def read(self,fields=None,selection=None):
        with open(self.path) as f:
            for line in f:
                line = line.strip().split()
                try:
                    chr,start,end = (line[0],int(line[1]),int(line[2]))
                except (IndexError,ValueError):
                    sys.exit(("Library 'bbcflib' not found. "
                              "Only 'bed' and 'bedGraph' formats available. "
                              "Wrong line in file %s:\n%s"
                              % (os.path.basename(self.path),'\t'.join(line)) ))
                yield tuple(line)
try:
    from bbcflib.btrack import track
    assert track
except ImportError:
    track = Parser

###############################################################################

class Memory(object):
    """Remembers what has already been read, so that one can go back."""
    def __init__(self):
        self.content = []

    def load(self,content):
        pass

    def get(self,content):
        pass

###############################################################################

class Reader(object):
    def __init__(self,trackList,nfeat,nbp,sel):
        self.tracks = [track(t) for t in trackList]
        self.sel = self.parse_selection(sel)
        self.chrom = self.init_chr()
        self.temp = []
        self.chrom_change = False
        self.ntimes = 1
        if nbp:
            self.nbp = nbp
            self.nfeat = None
        elif nfeat:
            self.nfeat = nfeat
            self.nbp = None

    def parse_selection(self,sel):
        if not sel: return None
        elif re.search('^chr[0-9XY]*:[0-9XY]+-[0-9XY]+',sel):
            chr,coord = sel.split(':')
            st,en = coord.split('-')
            return {'chr':chr,'start':(int(st),int(en)),'end':(int(st),int(en))}
        elif re.search('^chr[0-9XY]*$',sel):
            return {'chr':sel}
        else: print "Bad region formatting, got -s %s,." % sel; sys.exit(1)

    def init_chr(self):
        for t in self.tracks:
            try: return t.read().next()[0]
            except StopIteration: continue

    def read(self):
        """Yield a list of lists [[(1,2,n),(3,4,n)], [(1,3,n),(5,6,n)]] with either the *nfeat*
           next items, or all next items within an *nbp* window. `n` is a name or a score."""
        streams = [t.read(fields=t.fields[:4]) for t in self.tracks]
        self.go_to_selection(streams)
        if self.nfeat:
            content = self.read_nfeat(streams)
        elif self.nbp:
            content = self.read_nbp(streams)
        return content

    def go_to_selection(self,streams):
        nosel = (0,sys.maxint)
        skipped = 0
        self.temp = [s.next() for s in streams]
        if self.sel:
            for i,stream in enumerate(streams):
                try:
                    chrom,start,end = self.temp[i][:3]
                    while chrom != self.sel.get('chr',self.chrom) \
                    or not (self.sel.get('start',nosel)[0] < start < self.sel.get('start',nosel)[1]) \
                    or not (self.sel.get('end',nosel)[0] < end < self.sel.get('end',nosel)[1]):
                        self.temp[i] = stream.next()
                        chrom,start,end = self.temp[i][:3]
                        skipped += 1
                except StopIteration:
                    self.temp[i] = (self.chrom,0,0,'00')
            self.chrom = self.sel['chr']
            if self.nbp:
                temppos = [x[1] for x in self.temp if x[3]!='00']
                if temppos: self.ntimes += min(temppos) / self.nbp
                else: sys.exit("Chromosome %s not found." % self.chrom)
            elif self.nfeat:
                self.ntimes += skipped / self.nfeat

    def read_nfeat(self,streams):
        available_streams = range(len(streams))
        self.temp = [[x,k] for k,x in enumerate(self.temp)]
        sortkey = lambda x:(x[0][0],x[0][1])
        self.chrom_change = False
        # Repeat & yield each time the function is called
        while available_streams:
            self.chrom_change = False
            toyield = [[] for _ in streams]
            len_toyield = 0
            # Load *nfeat* in *toyield*
            while len_toyield < self.nfeat and self.temp:
                self.temp.sort(key=sortkey)
                min_idx = self.temp[0][-1]
                nextitem = self.temp.pop(0)[0]
                #if nextitem[2] >= minpos: # feature count is wrong
                #    toyield[min_idx].append(nextitem[1:])
                toyield[min_idx].append(nextitem[1:])
                len_toyield += 1
                try:
                    self.temp.append([streams[min_idx].next(),min_idx]) # read next item
                    chrom = self.temp[0][0][0]
                    if chrom != self.chrom:
                        self.chrom_change = True
                        self.chrom = chrom
                        break
                except StopIteration:
                    try: available_streams.pop(min_idx)
                    except IndexError: continue
            if any(toyield):
                # Add feats that go partially beyond
                if not self.chrom_change:
                    maxpos = max(x[-1][1] for x in toyield if x)
                    for n,x in enumerate(self.temp):
                        if x[0][1] < maxpos:
                            toyield[x[-1]].append((x[0][1],min(x[0][2],maxpos),x[0][3]))
                            if x[0][2] > maxpos:
                                self.temp[n][0] = (x[0][0],maxpos,x[0][2],x[0][3])
                #print "Toyield", toyield
                yield toyield
            else: break

    def read_nbp(self,streams):
        available_streams = range(len(streams))
        # Repeat & yield each time the function is called
        while available_streams:
            maxpos = self.ntimes*self.nbp
            toyield = [[] for _ in streams]
            toremove = []
            self.chrom_change = False
            chrom = [self.chrom for _ in streams]
            # Load items within *nbp* in *toyield*
            for n in available_streams:
                x = self.temp[n]
                while x[0] == self.chrom and x[2] <= maxpos:
                    toyield[n].append((x[1],x[2],x[3]))
                    try: x = streams[n].next()
                    except StopIteration:
                        toremove.append(n)
                        self.temp[n] = None
                        break
                if x[0] != self.chrom:
                    chrom[n] = x[0]
                elif x[2] > maxpos and x[1] < maxpos:
                    toyield[n].append((x[1],maxpos,x[3]))
                    x = (x[0],maxpos,x[2],x[3])
                self.temp[n] = x
            for n in toremove: available_streams.remove(n)
            if all(chrom[n] != self.chrom for n in available_streams):
                self.chrom = chrom[0]
                self.chrom_change = True
            if any(toyield):
                #print "Toyield", toyield
                yield toyield
            else:
                yield [[(0,0,'00')] for _ in streams]

###############################################################################

class Drawer(object):
    def __init__(self,names,types,nfeat,nbp):
        self.names = names # [file names]
        self.types = types # ['intervals' or 'density']
        self.nfeat = nfeat
        self.nbp = nbp
        self.ntimes = 0    # number of times the draw function is called
        self.maxpos = 0    # rightmost coordinate to display
        self.minpos = 0    # leftmost coordinate to display
        self.keydown = ''
        # Geometry
        self.root = tk.Tk()
        self.thisfeat = tk.Label() # popup showing the name of the feat under the mouse pointer
        self.WIDTH = 800   # window width
        self.htrack = 30   # canvas height
        self.rmargin = 100 # width of the right margin
        self.feat_pad = 10 # space between feat rectangles and the border of the canvas
        self.wlabel = 0    # width of the left margin with the track names
        self.wcanvas = 0   # width of the canvas
        self.reg_bp = 0    # size of the genomic region to display, in bp
        # Colors
        self.bg = "grey"
        self.canvas_bg = "white"
        self.feat_col = "blue"
        self.dens_col = "green"
        self.line_col = "black"

    def draw(self,content,chrom):
        """Create a new window and draw from the *content* coordinates
           (of the form [[(1,2,n),(3,4,n)], [(3,5,n),(4,6,n)],...],
           where `n` is either a name or a score)."""
        def keyboard(event):
            if event.keysym == 'Escape':
                self.keydown = chr(27)
                self.root.destroy()
                sys.exit(0)
            elif event.keysym == 'space':
                self.keydown = ' '
                self.root.quit()
            elif event.keysym == 'Left':
                self.keydown = chr(37)
                self.root.quit()
                pass # backwards??
            elif event.keysym == 'Right':
                self.keydown = chr(39)
                self.root.quit()
                pass # slowly forward
            elif event.keysym == 'BackSpace':
                self.keydown = chr(127)
                self.root.quit()
                pass # Return to the beginning
        self.root.title("gless")
        self.root.bind("<Key>", keyboard)
        self.root.config(bg=self.bg)
        self.minpos = self.maxpos if self.ntimes > 0 else 0
        if self.nbp:
            self.maxpos = self.ntimes*self.nbp
            self.reg_bp = self.maxpos - self.minpos
        elif self.nfeat:
            self.maxpos = max(t[-1][1] for t in content if t)
            self.reg_bp = self.maxpos - self.minpos
        self.reg_bp = float(max(self.reg_bp,self.nbp))
        self.draw_labels()
        self.draw_tracks(content)
        self.draw_margin(chrom)
        self.draw_axis(content)
        self.root.wm_attributes("-topmost", 1) # makes the window stay on top
        self.root.mainloop()

    def bp2px(self,x,wwidth,reg_bp):
        """Transform base pair coordinates to distances in pixels."""
        try: return x * wwidth/reg_bp
        except ZeroDivisionError: return 0

    def draw_labels(self):
        """Write track names on the left."""
        for n,name in enumerate(self.names):
            l = tk.Label(self.root,text=self.names[n],bd=0,highlightthickness=0,bg=self.bg,padx=5)
            self.wlabel = max(l.winfo_reqwidth(),self.wlabel)
            l.grid(row=n,column=0)
        self.wcanvas = self.WIDTH-self.wlabel-self.rmargin

    def draw_tracks(self,content):
        """Draw the canvas with the tracks in the middle."""
        def show_feat_name(event):
            canvas = event.widget
            x,y = event.x, event.y
            closest = canvas.find_closest(x,y)
            closest_coords = canvas.coords(closest)
            if closest_coords[0] != 0: # the base line has x=0 and is often closer
                self.thisfeat.place(x=closest_coords[0]+self.wlabel,
                                    y=closest_coords[1]+canvas.winfo_y(), anchor='w')
                self.thisfeat.config(text=name_map[canvas][closest[0]], bg='red')
                self.thisfeat.lift()
            else:
                self.thisfeat.place_forget()
        name_map = {} # correspondance canvas object id - feat name or score
        feat_thk = self.htrack - 2*self.feat_pad
        canvas = [tk.Canvas(self.root,height=self.htrack,bd=0,bg=self.canvas_bg,highlightthickness=0)
                  for _ in self.names]
        for n,t in enumerate(content):
            type = self.types[n]
            c = canvas[n]
            c.config(width=self.wcanvas)
            c.grid(row=n,column=1)
            c.bind("<Motion>", show_feat_name)
            name_map[c] = {}
            if type == 'intervals':
                c.create_line(0,self.htrack/2.,self.WIDTH,self.htrack/2.,fill=self.line_col) # track axis
                y1,y2 = (0+self.feat_pad,feat_thk+self.feat_pad)
                for k,feat in enumerate(t):
                    f1,f2,g = (feat[0],feat[1],feat[2])
                    x1 = self.bp2px(f1-self.minpos,self.wcanvas,self.reg_bp)
                    x2 = self.bp2px(f2-self.minpos,self.wcanvas,self.reg_bp)
                    if f1 == self.minpos: x1-=1 # no border
                    r = c.create_rectangle(x1,y1,x2,y2,fill=self.feat_col)
                    name_map[c][r] = g
            elif type == 'density':
                c.create_line(0,self.htrack-1,self.WIDTH,self.htrack-1,fill=self.line_col) # baseline
                if t:
                    top_bp = max(float(x[2]) for x in t) # highest score
                    top = self.bp2px(top_bp,self.htrack,top_bp)
                    for k,feat in enumerate(t):
                        f1,f2,s = (feat[0],feat[1],feat[2])
                        x1 = self.bp2px(f1-self.minpos,self.wcanvas,self.reg_bp)
                        x2 = self.bp2px(f2-self.minpos,self.wcanvas,self.reg_bp)
                        s = self.bp2px(float(s),self.htrack,top_bp)
                        if f1 == self.minpos: x1-=1 # no border
                        r = c.create_rectangle(x1,self.htrack-1,x2,self.htrack-s+5,fill=self.dens_col)
                        name_map[c][r] = str(s)
                    c.create_line(0,self.htrack-top+5,5,self.htrack-top+5) # vertical scale
                    c.create_line(2,self.htrack,2,self.htrack-top+5) # vertical scale
                    c.create_text(6,self.htrack-top+5,text=str(top_bp),anchor='w')

    def draw_margin(self,chrom):
        """Add a blank frame on the right as a margin, and the chromosome name."""
        w = tk.Label(text=chrom,bg='red')
        w.grid(row=0,column=2)
        for n in range(1,len(self.names)):
            w = tk.Frame(width=self.rmargin,height=self.htrack,bg=self.bg)
            w.grid(row=n,column=2)

    def draw_axis(self,content):
        """Draw the horizontal scale."""
        c = tk.Canvas(self.root,width=self.wcanvas,height=2*self.htrack,bd=0,
                      bg=self.canvas_bg,highlightthickness=0)
        c.grid(row=len(self.names)+1,column=1)
        pad = c.winfo_reqheight()/2.
        c.create_line(0,pad,self.WIDTH,pad,fill=self.line_col)  # axis
        for n,t in enumerate(content):
            for k,feat in enumerate(t):
                f1,f2 = (feat[0],feat[1])
                x1 = self.bp2px(f1-self.minpos,self.wcanvas,self.reg_bp)
                x2 = self.bp2px(f2-self.minpos,self.wcanvas,self.reg_bp)
                c.create_line(x1,pad,x1,pad-5,fill=self.line_col)
                c.create_line(x2,pad,x2,pad+5,fill=self.line_col)
                if f1!=self.minpos and f1!=self.maxpos:
                    c.create_text(x1,pad-5,text=str(f1),anchor='s')
                if f2!=self.minpos and f2!=self.maxpos:
                    c.create_text(x2,pad+5,text=str(f2),anchor='n')
        min_label = tk.Label(text=str(self.minpos),bd=0,bg=self.bg,anchor='e')
        min_label.grid(row=len(self.names)+1,column=0,sticky='e',padx=5)
        max_label = tk.Label(text=str(self.maxpos),bd=0,bg=self.bg,anchor='w')
        max_label.grid(row=len(self.names)+1,column=2,sticky='w',padx=5)

###############################################################################

class Gless(object):
    def __init__(self,trackList,nfeat,nbp,sel,fix):
        self.trackList = trackList
        self.names = [os.path.basename(t) for t in trackList]
        self.types = [self.get_type(n) for n in self.names]
        self.nfeat = nfeat
        self.nbp = nbp
        self.sel = sel
        self.fix = fix
        self.stream = None
        self.content = None
        self.needtodraw = True
        self.reader = Reader(self.trackList,self.nfeat,self.nbp,self.sel)
        self.drawer = Drawer(self.names,self.types,self.nfeat,self.nbp)

    def get_type(self,filename):
        """Return whether it is a track with 'intervals' or a 'density'."""
        with track(filename) as t:
            if t.format.lower() in ['bed']:
                return 'intervals'
            elif t.format.lower() in ['bedgraph','wig','bigWig','sga']:
                return 'density'

    def get_score_limits(self,fix):
        if not fix: return
        if len(fix.split(','))==1: return (0,float(fix))
        elif len(fix.split(','))==2: return tuple(fix.split(','))

    def __call__(self):
        """Main controller function."""
        self.stream = self.reader.read()
        try: self.content = self.stream.next()
        except StopIteration:
            print "Nothing to show"
            sys.exit(0)
        chrom = self.reader.chrom
        while True:
            if self.needtodraw:
                self.drawer.ntimes = self.reader.ntimes
                self.drawer.draw(self.content,chrom)
                self.needtodraw = False
                chrom = self.reader.chrom
            if self.drawer.keydown == chr(27): # "Esc" pressed: quit
                sys.exit(0)
            elif self.drawer.keydown == ' ': # "Space" pressed: next
                self.fast_forward()
            elif self.drawer.keydown == chr(127): # "BackSpace" ("Delete") pressed: return
                self.return_to_beginning()
                chrom = self.reader.chrom
            elif self.drawer.keydown == chr(37): # Left arrow pressed: shift left
                self.slow_reward()
            elif self.drawer.keydown == chr(39): # Right arrow pressed: shift right
                self.slow_forward()

    def reinit(self):
        self.drawer.minpos = 0
        self.drawer.maxpos = 0
        self.reader.ntimes = 1

    def load_next(self):
        try:
            self.content = self.stream.next() # Load next data
            for w in self.drawer.root.children.values(): # Clear the window
                w.destroy()
        except StopIteration:
            print "End of file"
            self.reader.ntimes -= 1
            self.drawer.ntimes -= 1
        self.needtodraw = True
        for w in self.drawer.root.children.values(): # Clear the window
            w.destroy()

    def return_to_beginning(self):
        self.reinit()
        self.reader = Reader(self.trackList,self.nfeat,self.nbp,self.sel)
        self.stream = self.reader.read()
        self.load_next()

    def fast_forward(self):
        if self.reader.chrom_change:
            self.reinit()
        else:
            self.reader.ntimes += 1
            self.drawer.ntimes += 1
        self.load_next()

    def slow_forward(self):
        self.load_next()

    def fast_reward(self):
        self.load_next()

    def slow_reward(self):
        self.load_next()

###############################################################################

"""If no explicit selection, chromosomes are base on the first track."""

def main():
    parser = argparse.ArgumentParser(description="Graphical 'less' for track files.")
    parser.add_argument('-n','--nfeat', default=10, type=int,
                       help="Number of features to display, exclusive with -b. [10]")
    parser.add_argument('-b','--nbp', default=None, type=int,
                       help="Number of base pairs to display, exclusive with -n.")
    parser.add_argument('-s','--sel', default=None,
                       help="Region to display, formatted as in 'chr1:12-34', or \
                             a chromosome only ('chr1').")
    parser.add_argument('-f','--fix', default=None,
                       help="Fixed range of scores for the vertical scale. One number \
                            (e.g. -f 10) indicates the max; two numbers separated by a comma \
                            (e.g. -f 5,10) indicate the min and max.")
    parser.add_argument('file', nargs='+', default=None,
                       help='A set of track files, separated by spaces')
    args = parser.parse_args()
    if args.nbp: args.nfeat = None
    elif args.nfeat and args.nfeat > 10000:
        print "Up to 10000 features permitted, got -n %s." % args.nfeat; sys.exit(1)
    Gless(args.file,args.nfeat,args.nbp,args.sel,args.fix)()

if __name__ == '__main__':
    sys.exit(main())

# python gless.py -b 20 testing_files/test1.bed testing_files/test2.bed testing_files/test1.bedgraph testing_files/test2.bedgraph testing_files/yeast_genes.bed


#fuckinfile = "/Users/julien/Library/Saved Application State/org.python.python.savedState"
#os.chmod(fuckinfile,0o777)
#os.remove(fuckinfile) # bug on osx

#trackList = ['testing_files/test1.bed','testing_files/test2.bed']
#root.focus_set()
#print c.winfo_reqheight(), c.winfo_reqwidth()
#print c.winfo_width(), c.winfo_height()
#bg = "#%02x%02x%02x" % (255, 255, 224) #beige background
    #c.create_line(0,0,0,self.htrack,fill="grey") # separator label|canvas
    #c.create_line(0,0,0,2*pad,fill=line_col)         # separator
        #if self.geometry: # keep previous position of the window
        #    self.geometry = '+'.join(['']+self.geometry.split('+')[1:]) # "+70+27"
        #    self.root.geometry(self.geometry)
