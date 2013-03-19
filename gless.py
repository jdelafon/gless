#!/usr/bin/env python
import Tkinter as tk
import os,sys
import argparse,re
import csv

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
            reader = csv.reader(f,delimiter='\t',quotechar='|')
            line0 = reader.next()
            nfields = len(line0)
            if not line0[0].startswith("track") or line0[0].startswith("#"): f.seek(0)
            fun = float if self.format.lower()=='bedgraph' else lambda x:x
            for line in reader:
                try:
                    chr,start,end = (line[0],int(line[1]),int(line[2]))
                    other = fun(line[3]) if nfields > 3 else "00"
                except (IndexError,ValueError):
                    sys.exit(("Library 'bbcflib' not found. "
                              "Only 'bed' and 'bedGraph' formats available. "
                              "Wrong line in file %s:\n%s"
                              % (os.path.basename(self.path),'\t'.join(line)) ))
                yield (chr,start,end,other)
try:
    from bbcflib.btrack import track
    assert track
except ImportError:
    print "Custom parser"
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
        self.available_streams = range(len(trackList))
        self.sel = sel
        self.temp = []
        self.chrom = self.init_chr()
        self.chrom_change = False
        self.next_chrom = self.chrom
        self.ntimes = 1
        if nbp:
            self.nbp = nbp
            self.nfeat = None
        elif nfeat:
            self.nfeat = nfeat
            self.nbp = None

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
        if self.sel and self.ntimes == 1:
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
                    self.temp.pop(i)
                    self.available_streams.remove(i)
            self.chrom = self.sel['chr']
            if self.nbp:
                temppos = [x[1] for x in self.temp if x[3]!='00']
                if temppos: self.ntimes += min(temppos) / self.nbp
                else: sys.exit("Chromosome %s not found." % self.chrom)
            elif self.nfeat:
                self.ntimes += skipped / self.nfeat

    def read_nfeat(self,streams):
        self.temp = [[x,n] for n,x in enumerate(self.temp)]
        toremove = []
        # Load *nfeat* feats of each track in the buffer
        for n in self.available_streams:
            for _ in range(self.nfeat):
                try:
                    self.temp.append([streams[n].next(),n])
                except StopIteration:
                    toremove.append(n)
                    break
        for n in toremove: self.available_streams.remove(n)
        # Repeat & yield each time the function is called
        while self.temp:
            self.chrom_change = False
            toyield = [[] for _ in streams]
            # Isolate one chromosome
            chrtemp = sorted([x for x in self.temp if x[0][0]==self.chrom], key=lambda x:x[0][2])
            rest = [x for x in self.temp if x[0][0]!=self.chrom]
            self.temp = chrtemp+rest
            if len(chrtemp) <= self.nfeat and rest:
                self.chrom_change = True
                self.next_chrom = rest[0][0][0]
            # Load *nfeat* in *toyield*
            for x,n in chrtemp[:self.nfeat]:
                toyield[n].append( x[1:3]+(x[3:] or ('00',)) )
                # Reload the buffer with one element for each element read,
                # so there are always *nfeat* x ntracks elements
                try: self.temp.append([streams[n].next(),n])
                except StopIteration:
                    try: self.available_streams.remove(n)
                    except ValueError: continue
            self.temp = self.temp[len(chrtemp[:self.nfeat]):]
            if any(toyield):
                # Add feats that go partially beyond
                maxpos = max(x[-1][1] for x in toyield if x)
                unseen = list(self.available_streams) # copy
                chrtemp = chrtemp[self.nfeat:]
                for k,feat in enumerate(chrtemp):
                    x,n = feat
                    if n in unseen and x[1] < maxpos:
                        toyield[n].append( (x[1],min(x[2],maxpos))+(x[3:] or ('00',)) )
                        self.temp[k][0] = (x[0],maxpos,x[2])+x[3:]
                        unseen.remove(n)
                    if not unseen: break
                #print "Toyield", toyield
                yield toyield
            else: break

    def read_nbp(self,streams):
        # Repeat & yield each time the function is called
        while self.available_streams:
            maxpos = self.ntimes*self.nbp
            toyield = [[] for _ in streams]
            toremove = []
            self.chrom_change = False
            chrom = [self.chrom for _ in streams]
            # Load items within *nbp* in *toyield*
            for n in self.available_streams:
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
            for n in toremove: self.available_streams.remove(n)
            if all(chrom[n] != self.chrom for n in self.available_streams):
                self.next_chrom = chrom[0]
                self.chrom_change = True
            if any(toyield):
                #print "Toyield", toyield
                yield toyield
            else:
                yield [[(0,0,'00')] for _ in streams]


###############################################################################

class Drawer(object):
    def __init__(self,names,types,nfeat,nbp,sel,fix):
        self.names = names # [file names]
        self.types = types # ['intervals' or 'density']
        self.nfeat = nfeat
        self.nbp = nbp
        self.sel = sel     # selection, of the type {'chr':'chr1','start':(1,1),'end':(2,2)}
        self.fix = fix     # (limits,) for the range of the vertical scale
        self.ntimes = 0    # number of times the draw function is called
        self.maxpos = 0    # rightmost coordinate to display
        self.minpos = 0    # leftmost coordinate to display
        self.keydown = ''
        # Geometry
        self.root = tk.Tk()
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
        self.root.focus_set() # not working?
        if self.nbp:
            if self.sel and self.sel.get('start') and self.maxpos == 0:
                self.minpos = self.sel['start'][0]
                self.maxpos = self.sel['end'][1]
            else:
                self.minpos = (self.ntimes-1)*self.nbp
                self.maxpos = self.ntimes*self.nbp
                self.reg_bp = self.maxpos - self.minpos
        elif self.nfeat:
            if self.ntimes > 1 and self.maxpos > 0:
                self.minpos = self.maxpos
                self.maxpos = max(t[-1][1] for t in content if t)
            elif self.sel and self.sel.get('start') and self.maxpos == 0:
                self.minpos = self.sel['start'][0]
                self.maxpos = self.sel['end'][1]
            else:
                self.minpos = max(0, min(t[0][0] for t in content if t))
                self.maxpos = max(t[-1][1] for t in content if t)
            self.reg_bp = self.maxpos - self.minpos
        self.reg_bp = float(max(self.reg_bp,self.nbp))
        self.draw_labels()
        self.draw_tracks(content)
        self.draw_rmargin(chrom)
        self.draw_axis(content)
        #self.root.wm_attributes("-topmost", 1) # makes the window stay on top
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
            if canvas.type(closest)=='rectangle':
                x1,y1,x2,y2 = canvas.coords(closest)
                # the base line has x=0 and is often closer
                self.thisfeat.place(x=x1+self.wlabel+(x2-x1)/2.,
                                    y=y1+canvas.winfo_y()+(y2-y1)/2., anchor='center')
                self.thisfeat.config(text=name_map[canvas][closest[0]],
                     bd=1,highlightbackground="black",highlightthickness=1)
                self.thisfeat.lift()
            else:
                self.thisfeat.place_forget()
        self.thisfeat = tk.Label() # popup showing the name of the feat under the mouse pointer
        name_map = {} # correspondance canvas object id - feat name or score
        feat_thk = self.htrack - 2*self.feat_pad
        canvas = [tk.Canvas(self.root,height=self.htrack,bd=0,bg=self.canvas_bg,highlightthickness=0)
                  for _ in self.names]
        for n,t in enumerate(content):
            type = self.types[n]
            c = canvas[n]
            c.config(width=self.wcanvas)
            c.grid(row=n,column=1,pady=5)
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
                hi = 2*self.htrack
                c.config(height=hi)
                c.create_line(0,hi-1,self.WIDTH,hi-1,fill=self.line_col) # baseline
                if t:
                    top_bp = self.fix[1] if self.fix else max(float(x[2]) for x in t) # highest score
                    top = self.bp2px(top_bp,hi,top_bp)
                    for k,feat in enumerate(t):
                        f1,f2,g = (feat[0],feat[1],feat[2])
                        x1 = self.bp2px(f1-self.minpos,self.wcanvas,self.reg_bp)
                        x2 = self.bp2px(f2-self.minpos,self.wcanvas,self.reg_bp)
                        s = float(g)-self.fix[0] if self.fix else float(g)
                        s = self.bp2px(s,hi,top_bp)
                        if f1 == self.minpos: x1-=1 # no border
                        if s > 0:
                            r = c.create_rectangle(x1,hi-1,x2,hi-s+5,fill=self.dens_col)
                            name_map[c][r] = str(g)
                    c.create_line(0,hi-top+5,5,hi-top+5) # vertical scale
                    c.create_line(2,hi,2,hi-top+5) # vertical scale
                    c.create_text(6,hi-top+5,text=str(top_bp),anchor='w')
        back = tk.Frame(self.root,bg=self.canvas_bg,width=self.wcanvas) # blank background
        back.grid(column=1,row=0,rowspan=len(self.names),sticky=["N","S"])
        back.lower()

    def draw_rmargin(self,chrom):
        """Add a blank frame on the right as a margin, and the chromosome name."""
        w = tk.Label(text=chrom,bg='white')
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
                # Need to not write overlapping pos
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
        self.nfeat = nfeat
        self.nbp = nbp
        self.sel = self.parse_selection(sel)
        self.names = [os.path.basename(t) for t in trackList]
        self.types = [self.get_type(n) for n in self.names]
        self.stream = None
        self.content = None
        self.needtodraw = True
        fix = self.get_score_limits(fix)
        self.reader = Reader(self.trackList,self.nfeat,self.nbp,self.sel)
        self.drawer = Drawer(self.names,self.types,self.nfeat,self.nbp,self.reader.sel,fix)

    def get_type(self,filename):
        """Return whether it is a track with 'intervals' or a 'density'."""
        with track(filename) as t:
            if t.format.lower() in ['bed']:
                return 'intervals'
            elif t.format.lower() in ['bedgraph','wig','bigWig','sga']:
                return 'density'

    def parse_selection(self,sel):
        if not sel: return None
        elif re.search('^chr[0-9XY]*:[0-9XY]+-[0-9XY]+',sel):
            chr,coord = sel.split(':')
            st,en = coord.split('-')
            return {'chr':chr,'start':(int(st),int(en)),'end':(int(st),int(en))}
        elif re.search('^chr[0-9XY]*$',sel):
            return {'chr':sel}
        else: print "Bad region formatting, got -s %s,." % sel; sys.exit(1)

    def get_score_limits(self,fix):
        if not fix: return
        elif len(fix.split(','))==1: return (0,float(fix))
        elif len(fix.split(','))==2: return tuple(float(x) for x in fix.split(','))
        else: sys.exit("Wrong format for -f option: got %s." % fix)

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
                if self.reader.chrom_change:
                    self.reader.chrom = self.reader.next_chrom
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
            if self.nfeat:
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



#fuckinfile = "/Users/julien/Library/Saved Application State/org.python.python.savedState"
#os.chmod(fuckinfile,0o777)
#os.remove(fuckinfile) # bug on osx

#print c.winfo_reqheight(), c.winfo_reqwidth()
#print c.winfo_width(), c.winfo_height()
#bg = "#%02x%02x%02x" % (255, 255, 224) #beige background
    #c.create_line(0,0,0,self.htrack,fill="grey") # separator label|canvas
    #c.create_line(0,0,0,2*pad,fill=line_col)         # separator
        #if self.geometry: # keep previous position of the window
        #    self.geometry = '+'.join(['']+self.geometry.split('+')[1:]) # "+70+27"
        #    self.root.geometry(self.geometry)
