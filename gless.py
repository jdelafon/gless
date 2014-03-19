#!/usr/bin/env python

"""
Description
===========
This application is designed to be used as a graphical equivalent of the `less` command
in Unix to better visualize track files. It will give you an insight of the content of
your files that is probably telling more than columns of numbers.

If the library `bbcflib` is found on your system, `track` will be used and will
recognize .bed, .bedgraph, .wig, .sga, .bigWig, .sql, .sam formats. Else it can still
read .bed and .bedGraph files.

Since files are read sequentially (without loading temp in memory), it does not read
backwards. Also the chromosomes names present in each file are not known in advance.
The first file you give as input is taken as a reference. This can cause chromosomes
to be skipped in the secondary files, if their names or order is different.

Usage
=====
Press the SPACE bar to read forward, RETURN (or Delete) to return to the beginning,
ESC to quit. Move the cursor on elements to display their name/score.

Options:

* -n nfeat: display the next *nfeat* features (from all tracks together).
* -b nbp: display the next *nbp* base pairs window.
* -y ylim: set the vertical scale for numeric tracks: either `<min>,<max>` or just `<max>`.
* -s sel: selection: either a chromosome name, or a region specified as <chr>:<start>.
  The right bound is set by the -n/-b argument."

Known issues:
=============

* Bug on OSX: if it says something like "Could not restore the previous window", remove
  `/Users/<User>/Library/Saved Application State/org.python.python.savedState` .
* Focus on the main window is not automatic on some OX.
* The main window may not appear on top on some OS.
"""

import Tkinter as tk
import os,sys
import argparse,re
import csv

###############################################################################

class Parser(object):
    #A replacement for track when bbcflib is not found, able to parse
    #bed and bedGraph formats only. Called in the Reader class."""
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
    from bbcflib.track import track
    assert track
except ImportError:
    track = Parser

###############################################################################

class Memory(object): # Not working yet
    def __init__(self):
        self.content = []
        self.chrom = None
        self.limit = 10

    def save(self,content):
        chrom = content[0][0][0]
        self.content.setdefault(chrom,[]).append(content)

    def load(self):
        return self.content[-1]

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
        """Find the initial chromosome name."""
        for t in self.tracks:
            try: return t.read().next()[t.fields.index('chr')]
            except StopIteration: continue

    def read(self):
        """Yield a list of lists [[(1,2,n),(3,4,n)], [(1,3,n),(5,6,n)]] with either the *self.nfeat*
           next items, or all next items within an *self.nbp* window. `n` is a name or a score."""
        streams = []
        for t in self.tracks:
            if all(f in t.fields for f in ["chr","start","end"]):
                _f = ["chr","start","end"]
                if "name" in t.fields:
                    _f.append("name")
            else:
                _f = t.fields[:4]
            streams.append(t.read(fields=_f))
        self.go_to_selection(streams)
        if self.nfeat:
            content = self.read_nfeat(streams)
        elif self.nbp:
            content = self.read_nbp(streams)
        return content

    def go_to_selection(self,streams):
        """Skip all features not passing the selection filter before filling the buffer."""
        skipped = 0
        self.temp = [s.next() for s in streams]
        if self.sel and self.ntimes == 1:
            selected_chrom = self.sel.get('chr',self.chrom)
            selected_start = self.sel.get('start',[0])[0]
            for i,stream in enumerate(streams):
                try:
                    chrom,start,end = self.temp[i][:3]
                    while chrom != selected_chrom:
                        self.temp[i] = stream.next()
                        chrom,start,end = self.temp[i][:3]
                        skipped += 1
                    while start < selected_start:
                        if end > selected_start: # overlapping
                            self.temp[i] = (chrom,selected_start,end)+self.temp[i][3:]
                            break
                        else:
                            self.temp[i] = stream.next()
                            chrom,start,end = self.temp[i][:3]
                            skipped += 1
                except StopIteration:
                    self.temp.pop(i)
                    self.available_streams.remove(i)
                except IndexError:
                    sys.exit("Unknown region.")
            self.chrom = self.sel['chr']
            if self.nbp:
                temppos = [x[1] for x in self.temp if x[3]!='00']
                if temppos: self.ntimes += min(temppos) / self.nbp
                else: sys.exit("Chromosome %s not found." % self.chrom)
            elif self.nfeat:
                self.ntimes += skipped / self.nfeat

    def read_nfeat(self,streams):
        """Yield the next *nfeat* features."""
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
                yield toyield
            else: break

    def read_nbp(self,streams):
        """Yield all features in the next *nbp* base pairs window."""
        # Repeat & yield each time the function is called
        shift = self.sel.get('start',[0])[0] if self.sel else 0
        while self.available_streams:
            maxpos = self.ntimes*self.nbp + shift
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
                yield toyield
            else:
                yield [[(0,0,'00')] for _ in streams]

###############################################################################

class Drawer(object):
    def __init__(self,names,types,nfeat,nbp,sel,ylim):
        self.names = names # [file names]
        self.types = types # ['intervals' or 'density']
        self.nfeat = nfeat
        self.nbp = nbp
        self.sel = sel     # selection, of the type {'chr':'chr1','start':(1,1),'end':(2,2)}
        self.ylim = ylim   # (limits,) for the range of the vertical scale
        self.ntimes = 0    # number of times the draw function is called
        self.maxpos = 0    # rightmost coordinate to display
        self.minpos = 0    # leftmost coordinate to display
        self.nticks = 10   # number of ticks on the horiz axis if regular scale
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

        def set_boundaries():
            if self.nbp:
                # first time: shift to selection
                if self.sel and self.sel.get('start') and self.maxpos == 0:
                    self.minpos = self.sel['start'][0]
                    self.maxpos = self.minpos + self.nbp
                # next times: keep this shift
                elif self.sel and self.sel.get('start'):
                    self.minpos = (self.ntimes-1)*self.nbp + self.sel['start'][0]
                    self.maxpos = self.ntimes*self.nbp + self.sel['start'][0]
                # no selection
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
                else:
                    self.minpos = max(0, min(t[0][0] for t in content if t))
                self.maxpos = max(t[-1][1] for t in content if t)
                self.reg_bp = self.maxpos - self.minpos
            self.reg_bp = float(max(self.reg_bp,self.nbp))

        self.root.title("gless")
        self.root.bind("<Key>", keyboard)
        self.root.config(bg=self.bg)
        self.root.focus_set() # not working?
        self.root.resizable(0,0) # disable window resizing
        set_boundaries()
        self.draw_labels()
        self.draw_tracks(content)
        self.draw_rmargin(chrom)
        self.draw_axis(content)
        try: self.root.wm_attributes("-topmost", 1) # makes the window stay on top
        except: pass # depends on the OS
        def _finish():
            self.root.destroy()
            sys.exit(0)
        self.root.protocol("WM_DELETE_WINDOW", _finish)
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
                    if g == '00': g = "%d-%d" % (f1,f2)
                    name_map[c][r] = g
            elif type == 'density':
                m = 6 # margin
                hi = 2 * self.htrack # twice higher than for intervals
                c.config(height=hi)
                if t:
                    ymax = self.ylim.get('max', max(float(x[2]) for x in t)) # max score to display
                    ymin = self.ylim.get('min', min(float(x[2]) for x in t)) # min score to display
                    yrange = abs(ymax-ymin)
                    if yrange: scale = (hi-2*m) / yrange # score to px
                    else: scale = 1
                    mid = max(0, ymax*scale +m-1)
                    for k,feat in enumerate(t):
                        f1,f2,g = (feat[0],feat[1],feat[2])
                        x1 = self.bp2px(f1-self.minpos,self.wcanvas,self.reg_bp)
                        x2 = self.bp2px(f2-self.minpos,self.wcanvas,self.reg_bp)
                        s = float(g)
                        if f1 == self.minpos: x1-=1 # no border
                        if s > 0:
                            if s > self.ylim.get('min',-1) > 0:
                                s = s - self.ylim['min']
                            if self.ylim.get('max'):
                                s = min(s,self.ylim['max'])
                            spx = mid-s*scale -1
                        else:
                            if s < self.ylim.get('max',1) < 0:
                                s = s - self.ylim['max']
                            if self.ylim.get('min'):
                                s = max(s,self.ylim['min'])
                            spx = mid-s*scale +1
                        r = c.create_rectangle(x1,mid,x2,spx,fill=self.dens_col)
                        name_map[c][r] = str(g)
                    ymax_px = m
                    ymin_px = yrange*scale + m
                    if ymax > 0:
                        c.create_line(0,ymax_px,5,ymax_px) # little horizontal tick, max
                        c.create_text(6,ymax_px,text=str(ymax),anchor='w') # max label
                    if ymin < 0:
                        c.create_line(0,ymin_px,5,ymin_px) # little horizontal tick, min
                        c.create_text(6,ymin_px,text=str(ymin),anchor='w') # min label
                    c.create_line(2,max(ymin_px,mid),2,min(mid,ymax_px)) # vertical scale
                    bl = min(hi-1,max(0, ymax*scale +m-1)) # position of the baseline
                    c.create_line(0,bl,self.wcanvas,bl,fill=self.line_col) # baseline
                else:
                    c.create_line(0,hi/2,self.wcanvas,hi/2,fill=self.line_col,dash=1) # baseline
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
        if sum(len(c) for c in content) <= 10*len(content):
            for n,t in enumerate(content):
                for k,feat in enumerate(t):
                    f1,f2 = (feat[0],feat[1])
                    x1 = self.bp2px(f1-self.minpos,self.wcanvas,self.reg_bp)
                    x2 = self.bp2px(f2-self.minpos,self.wcanvas,self.reg_bp)
                    # ticks
                    c.create_line(x1,pad,x1,pad-5,fill=self.line_col)
                    c.create_line(x2,pad,x2,pad+5,fill=self.line_col)
                    # labels
                    if f1!=self.minpos and f1!=self.maxpos:
                        c.create_text(x1,pad-5,text=str(f1),anchor='s')
                    if f2!=self.minpos and f2!=self.maxpos:
                        c.create_text(x2,pad+5,text=str(f2),anchor='n')
        else: # regular, linear scale
            ticksize = (self.maxpos-self.minpos) // self.nticks or 1
            for n,k in enumerate(range(self.minpos+ticksize,self.maxpos,ticksize)):
                x = self.bp2px(k-self.minpos,self.wcanvas,self.reg_bp)
                if n%2 == 0:
                    c.create_line(x,pad,x,pad-5,fill=self.line_col)
                    c.create_text(x,pad-5,text=str(k),anchor='s')
                else:
                    c.create_line(x,pad,x,pad+5,fill=self.line_col)
                    c.create_text(x,pad+5,text=str(k),anchor='n')
        min_label = tk.Label(text=str(self.minpos),bd=0,bg=self.bg,anchor='e')
        min_label.grid(row=len(self.names)+1,column=0,sticky='e',padx=5)
        max_label = tk.Label(text=str(self.maxpos),bd=0,bg=self.bg,anchor='w')
        max_label.grid(row=len(self.names)+1,column=2,sticky='w',padx=5)

###############################################################################

class Gless(object):
    def __init__(self,trackList,nfeat,nbp,sel,ylim):
        self.trackList = trackList
        self.nfeat = nfeat
        self.nbp = nbp
        self.sel = self.parse_selection(sel)
        self.names = [os.path.basename(t) for t in trackList]
        self.types = [self.get_type(n) for n in self.names]
        self.stream = None
        self.content = None
        self.needtodraw = True
        ylim = self.get_score_limits(ylim)
        self.reader = Reader(self.trackList,self.nfeat,self.nbp,self.sel)
        self.drawer = Drawer(self.names,self.types,self.nfeat,self.nbp,self.reader.sel,ylim)
        self.memory = Memory() # Not working yet

    def get_type(self,filename):
        """Return whether it is a track with 'intervals' or a 'density'."""
        with track(filename) as t:
            if t.format.lower() in ['bed','sam','bam']:
                return 'intervals'
            elif t.format.lower() in ['bedgraph','wig','bigWig','sga']:
                return 'density'

    def parse_selection(self,sel):
        """Transform 'chr1:12' into {'chr':'chr1','start':(12,12)}."""
        if not sel: return None
        elif re.search('^chr[0-9XY]*:[0-9XY]+',sel):
            chr,start = sel.split(':')
            return {'chr':chr,'start':(int(start),int(start))}
        elif re.search('^chr[0-9XY]*$',sel):
            return {'chr':sel}
        else: sys.exit("Bad region formatting, got '-s %s' ." % sel)

    def get_score_limits(self,ylim):
        """Transform '5,100' into {'min':5,'max':100}."""
        if not ylim: return {}
        elif len(ylim.split(','))==1: return {'max':float(ylim)}
        elif len(ylim.split(','))==2: return {'min':float(ylim.split(',')[0]),
                                              'max':float(ylim.split(',')[1])}
        else: sys.exit("Wrong format for -f option: got %s." % ylim)

    def __call__(self):
        """Main controller function."""
        self.stream = self.reader.read()
        try: self.content = self.stream.next()
        except StopIteration:
            sys.exit("Nothing to show")
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
        """Called after chrom change or returning to the beginning."""
        self.drawer.minpos = 0
        self.drawer.maxpos = 0
        self.reader.ntimes = 1

    def load_next(self):
        """Load next set of features and draw the new figure."""
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

def main():
    parser = argparse.ArgumentParser(description="Graphical 'less' for track files\n. \
                       Press the SPACE bar to read forward, RETURN (or Delete) to \
                       return to the beginning, ESC to quit.")
    parser.add_argument('-n','--nfeat', default=10, type=int,
                       help="Number of features to display, exclusive with -b. [10]")
    parser.add_argument('-b','--nbp', default=None, type=int,
                       help="Number of base pairs to display, exclusive with -n.")
    parser.add_argument('-s','--sel', default=None,
                       help="Region to display, formatted as <chr>:<start> (e.g. 'chr1:12'),\
                             or a chromosome name only ('chr1'). The right bound is set \
                             by the -n/-b argument.")
    parser.add_argument('-y','--ylim', default=None,
                       help="Fixed range of scores for the vertical scale. One number \
                            (e.g. -y 10) indicates the max positive value to display; \
                            two numbers separated by a comma indicate the min and the max. \
                            For negative values, make sure to use the equal sign (e.g. -y=-5,10).")
    parser.add_argument('file', nargs='+', default=None,
                       help='A set of track files, separated by spaces')
    args = parser.parse_args()
    if args.nbp: args.nfeat = None
    Gless(args.file,args.nfeat,args.nbp,args.sel,args.ylim)()

if __name__ == '__main__':
    sys.exit(main())

#------------------------------------------------------#
# This code was written by Julien Delafontaine         #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# julien.delafontaine@yandex.com                       #
#------------------------------------------------------#
