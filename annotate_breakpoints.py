#!/usr/bin/python

from Tkinter import *
import os,csv

## profile http://docs.python.org/library/profile.html
## import cProfile
## cProfile.run('foo()', 'fooprof')
## run snake run
## http://docs.python.org/library/multiprocessing.html
class Profiles(dict):
    def __init__(self,fn):
        self.read(fn)
    def get_profile(self,profile_id):
        if profile_id not in self:
            print "Loading profile", profile_id # to show status while loading
            self[profile_id]=Profile(profile_id)
        return self[profile_id]
    def read(self,fn):
        """Read csv profile table."""
        f=open(fn)
        #lines=f.readlines(3000000)
        reader = csv.reader(f,lineterminator="\n",delimiter=",")
        colnames = reader.next()
        colnums = dict([(n,i) for i,n in enumerate(colnames)])
        cid_col = colnums["profile.id"]
        for row in reader:
            ##d=dict(zip(colnames,row))
            p=self.get_profile(row[cid_col])
            p.add_point(*row)
        f.close()
        self.CHROMOSOME_ORDER = [str(x) for x in range(1,23)+["X"]]
        # all points are read, now cache values
        for p in self.values():
            p.cache_minmax()
        ## store min and max positions for each chrom
        for fun_name in "min","max":
            k="position"+fun_name
            fun = eval(fun_name)
            vals=[]
            for cname in self.CHROMOSOME_ORDER:
                items=[p[cname][k] for p in self.values() if cname in p]
                if items:
                    vals.append(fun(items))
            setattr(self,k,vals)
        self.chrom_lengths = [
            l-f for f,l in zip(self.positionmin,self.positionmax)]
        ## total length of all chromosomes
        self.total_chrom_length = sum([float(x) for x in self.chrom_lengths])
        #print self.chrom_lengths,self.CHROMOSOME_ORDER
        ## proportion of each chromosome
        self.chrom_lprops = [
            x/self.total_chrom_length for x in self.chrom_lengths]

class Profile(dict):
    def get_importance(self):
        return self.importance
    def __init__(self,profile_id):
        self.profile_id = profile_id
        self.chrom_order = []
    def add_point(self,profile_id,chromosome,position,logratio,*args,**kwargs):
        if chromosome not in self:
            self.chrom_order.append(chromosome)
            self[chromosome]={
                "position":[],
                "logratio":[],
                }
        c=self[chromosome]
        c["position"].append(int(position))
        c["logratio"].append(float(logratio))
    def cache_minmax(self):
        for c in self.values(): ## each chrom
            for k in 'position','logratio':
                for fun_name in "min","max":
                    fun = eval(fun_name)
                    c[k+fun_name]=fun(c[k])
            c["n"]=len(c[k])
        for fun_name in "min","max":
            k="logratio"+fun_name
            fun = eval(fun_name)
            val=fun([c[k] for c in self.values()])
            setattr(self,k,val)
        self.logratio_range = self.logratiomax - self.logratiomin
    def __repr__(self):
        items=[str(self.get(k,{"n":"."})["n"]) for k in self.chrom_order]
        return "Profile: "+self.profile_id+" ".join(items)

class Region(dict):
    def __init__(self,annotation,first,last):
        self.update({
                "annotation":annotation,
                "first":first,
                "last":last,
                "size":last-first,
                })
    def onClick(self,e):
        ## cycle through annotations as specified in ANNOTATION_ORDER
        if self["annotation"] in NEXT_ANNOTATION:
            self["annotation"]=NEXT_ANNOTATION[self["annotation"]]
            fill=REGION_COLORS[self["annotation"]]
            e.widget.itemconfig(self.id,fill=fill,outline=fill)
        else:
            e.widget.delete(self.id)
            self.regions.pop(self.region_index)
        return "break"

## this defines the colors to display for the annotations and the
## order in which to display them
ANNOTATION_ORDER=[
        ("1breakpoint","#ff7d7d"),
        ("0breakpoints",'#f6f4bf'),
        ]
ANNOTATIONS=[a for a,c in ANNOTATION_ORDER]
NEXT_ANNOTATION=dict(zip(ANNOTATIONS,ANNOTATIONS[1:]))
REGION_COLORS=dict(ANNOTATION_ORDER)

class RegionList(dict):
    """container for regions.

    facilitates deletion more automatically by passing the keys and
    the container reference to the items themselves, so they can
    remove themselves."""
    def __init__(self):
        self.counter = 0
    def add(self,r):
        r.regions = self
        self[self.counter] = r
        r.region_index = self.counter
        self.counter += 1

class ChromDisplay(Canvas):
    """Canvas with special onclick"""
    def onClick(self,e):
        self.orig_x = self.canvasx(e.x)
        self.resize_rect(e)
    def onMotion(self,e):
        self.delete(self.new_id)
        self.resize_rect(e)
    def resize_rect(self,e):
        x = self.canvasx(e.x)
        if x < 0:
            x=0
        if x > self.w:
            x=self.w
        if self.orig_x < x:
            self.left = self.orig_x
            self.right = x
        else:
            self.left = x
            self.right = self.orig_x
        self.new_id = self.make_rect(
            self.left,self.right,ANNOTATION_ORDER[0][1])
    def make_rect(self,left,right,fill):
        id = self.create_rectangle(
            left,1,right,self.h,fill=fill,outline=fill,activeoutline="black",
            tag="region")
        self.tag_lower("lines")
        self.tag_lower("region")
        self.tag_lower(self.bgid)
        self.tag_lower("interval")
        return id
    def to_position(self,pixels):
        return int(pixels * self.l / self.w + self.m)
    def to_pixels(self,position):
        #print position,self.m,self.w,self.l
        return int(float(position-self.m)*self.w/self.l)
    def onRelease(self,e):
        r = Region(
            ANNOTATION_ORDER[0][0],
            self.to_position(self.left),
            self.to_position(self.right),
            )
        r.id = self.new_id
        self.tag_bind(self.new_id,"<Button-1>",r.onClick)
        self.regions.add(r)
        #print self.regions
    def doneAnnotating(self,e):
        """Remove this profile and add another to the display."""
        ann = self.annotator
        importance = [
            (p.get_importance(),k) for k,p in ann.profiles.iteritems()]
        importance.sort(key=lambda t: t[0],reverse=True)
        #print importance
        p = ann.profiles[importance[0][1]]
        ann.bind_profile(p,self.row)

class Annotator(object):
    def onClose(self):
        """Save annotations to file before quitting."""
        lines = ["profile.id,chromosome,min,max,annotation"]
        for j,chr in enumerate(self.profiles.CHROMOSOME_ORDER):
            m = self.profiles.positionmin[j]
            M = self.profiles.positionmax[j]
            for pid in self.all_pids:
                region_list = self.regions_dict[(pid,chr)].values()
                for r in region_list:
                    lines.append("%s,%s,%d,%d,%s"%(
                            pid,
                            chr,
                            r["first"],
                            r["last"],
                            r["annotation"],
                            ))
        text = "\n".join(lines)
        f=open(self.regions_file,"w")
        f.write(text)
        f.close()
        self.root.destroy()
    def __init__(self,root,profiles_file,regions_file):
        self.root = root
        root.protocol("WM_DELETE_WINDOW",self.onClose)
        ## save file names for later to save annotations
        self.profiles = Profiles(profiles_file)
        self.all_pids = [p.profile_id for p in self.profiles.values()]
        self.active_profiles = {}

        self.regions_dict={}
        for chr in self.profiles.CHROMOSOME_ORDER:
            for pid in self.profiles:
                self.regions_dict[(pid,chr)] = RegionList()
        self.regions_file = regions_file
        ## if file does not exist, will create when we exit
        if os.path.isfile(regions_file): 
            f=open(regions_file)
            reader=csv.reader(f,lineterminator="\n",delimiter=",")
            reader.next()
            for p,c,first,last,annotation in reader:
                r = Region(annotation,int(first),int(last))
                k = (p,c)
                if k not in self.regions_dict:
                    self.regions_dict[k] = RegionList()
                self.regions_dict[k].add(r)
            f.close()

        MIN_CHROM_HEIGHT = 90
        WIDTH=root.winfo_screenwidth()-100
        HEIGHT=root.winfo_screenheight() -100
        ROWS = HEIGHT / MIN_CHROM_HEIGHT
        CHROM_HEIGHT = HEIGHT / ROWS - 2
        #print WIDTH,HEIGHT,ROWS,CHROM_HEIGHT,self.profiles.CHROMOSOME_ORDER
        root.geometry("%sx%s+0+0"%(WIDTH,HEIGHT))
        self.canvases = {}
        for j,c in enumerate(self.profiles.CHROMOSOME_ORDER):
            w = int(self.profiles.chrom_lprops[j] * WIDTH)-2
            for i in range(ROWS):
                widget = ChromDisplay(root,background="white",width=w,
                                height=CHROM_HEIGHT,borderwidth=0,
                                highlightthickness=1,
                                      highlightbackground="grey",
                                      highlightcolor="red",
                                      )
                ## all this stuff doesnt change with profiles (but
                ## will change with resizing)
                widget.m = self.profiles.positionmin[j]
                widget.M = self.profiles.positionmax[j]
                widget.l = self.profiles.chrom_lengths[j]
                widget.w = w
                widget.h = CHROM_HEIGHT
                #print widget.m,widget.M,widget.l,widget.w,widget.h
                widget.annotator = self ## for right-click callbacks
                widget.row = i
                self.canvases[(i,c)] = widget
                widget.grid(row=i,column=j,padx=0,pady=0)
        self.active_pids = range(ROWS)
        self.bound_time = 0
        #self.model = Model()
        for i,p in enumerate(self.profiles.values()):
            p.bind_count = 0
            p.last_time_shown = 1
            p.importance = 1.0
            if i < ROWS:
                self.bind_profile(p,i)
    def bind_profile(self,p,i):
        """Bind profile p to display row i."""
        ## first save profile that was here
        pid_earlier = self.active_pids[i]
        if pid_earlier in self.active_profiles:
            profile_earlier = self.active_profiles[pid_earlier]
            profile_earlier.last_time_shown = self.bound_time
            profile_earlier.importance = 1.0/self.bound_time
            self.profiles[pid_earlier] = profile_earlier
        ## then add the specified profile
        pid = p.profile_id
        p.bind_count += 1
        self.bound_time += 1
        p.bound_time = self.bound_time
        self.profiles.pop(pid)
        self.active_profiles[pid] = p
        self.active_pids[i] = pid
        for c in self.profiles.CHROMOSOME_ORDER:
            k=(pid,c)
            w=self.canvases[(i,c)]
            w.bind("<Button-3>",w.doneAnnotating)
            w.delete(ALL) ## delete previous contents first
            LINES=[0,-1,1]
            LINES_PX = [
                w.h * (1 -(x-p.logratiomin)/p.logratio_range)
                for x in LINES
                ]
            y=LINES_PX[0]
            w.create_rectangle(0,LINES_PX[1],w.w+1,LINES_PX[2],
                               fill="#e5e5e5",outline="",
                               tag="interval")
            w.bgid = w.create_rectangle(0,0,w.w,w.h,outline="",fill="")
            w.create_line(0,y,w.w,y,tag="lines")
            w.tag_bind(w.bgid,"<Button-1>",w.onClick)
            w.tag_bind(w.bgid,"<B1-Motion>",w.onMotion)
            w.tag_bind(w.bgid,"<ButtonRelease-1>",w.onRelease)
            w.regions = self.regions_dict[k]
            for r in w.regions.values():
                first = w.to_pixels(r["first"])
                last = w.to_pixels(r["last"])
                fill = REGION_COLORS[r["annotation"]]
                r.id = w.make_rect(first,last,fill)
                w.tag_bind(r.id,"<Button-1>",r.onClick)
            if c in p: ## for entire unobserved chromosomes
                chrom = p[c]
                norm = [
                    (lr-p.logratiomin)/p.logratio_range
                    for lr in chrom["logratio"]
                    ]
                y_px = [
                    w.h - n*w.h
                    for n in norm
                    ]
                x_px = [
                    (pos-float(w.m))*w.w/w.l
                    for pos in chrom["position"]
                    ]
                for x,y in zip(x_px,y_px):
                    w.create_oval(x-1,y-1,x+1,y+1,fill="",outline="blue")

if __name__ == "__main__":
    import sys
    root = Tk()
    #root.state("zoomed") # to start maximized
    args=sys.argv[1:]
    if len(args) != 2:
        print """Usage: %s profiles.csv annotations.csv
profiles.csv has columns profile.id,chromosome,position,logratio
annotations.csv has columns profile.id,chromosome,min,max,annotation
"""%sys.argv[0]
        sys.exit(1)
    ann = Annotator(root,*args)
    root.mainloop()
