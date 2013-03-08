import Tkinter as tk

def graph_points(seq, width=375, height=325):
    root = tk.Tk()
    c = tk.Canvas(root, width=width, height=height, bg='white')
    c.pack()
    y_stretch = 15
    y_gap = 20
    x_stretch = 10
    x_width = 20
    x_gap = 20
    for x, y in enumerate(data):
        x0 = x * x_stretch + x * x_width + x_gap
        y0 = height - (y * y_stretch + y_gap)
        x1 = x * x_stretch + x * x_width + x_width + x_gap
        y1 = height - y_gap
        c.create_rectangle(x0, y0, x1, y1, fill="red")
        c.create_text(x0+2, y0, anchor=tk.SW, text=str(y))
    root.mainloop()

data = (18, 15, 10, 7, 5, 4, 2, 5, 8, 10, 13)
graph_points(data)

