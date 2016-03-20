import numpy as np

# Credit to Steven Lambert for the tutorial

class Rectangle:

    def __init__(self, x: float, y: float, width: float, height: float):
        self.x = float(x)
        self.y = float(y)
        self.width = float(width)
        self.height = float(height)


class Quadtree:

    max_objects = 10
    max_levels = 5

    def __init__(self, p_level: int, p_bounds: Rectangle):
        self.level = int(p_level)
        self.objects = []
        self.bounds = p_bounds
        self.nodes = [None, None, None, None]

    def clear(self):
        self.objects.clear()
        for node in self.nodes:
            if node is not None:
                node.clear()

    def split(self):
        subwidth = self.bounds.width/2
        subheight = self.bounds.height/2
        x = self.bounds.x
        y = self.bounds.y

        self.nodes[0] = Quadtree(self.level + 1, Rectangle(x + subwidth, y, subwidth, subheight))
        self.nodes[1] = Quadtree(self.level + 1, Rectangle(x, y, subwidth, subheight))
        self.nodes[2] = Quadtree(self.level + 1, Rectangle(x, y + subheight, subwidth, subheight))
        self.nodes[3] = Quadtree(self.level + 1, Rectangle(x + subwidth, y + subheight, subwidth, subheight))

    def get_index(self, p_rect: Rectangle):
        index = -1
        v_midpoint = self.bounds.x + self.bounds.width/2
        h_midpoint = self.bounds.y + self.bounds.height/2

        top_quad = p_rect.y < h_midpoint and p_rect.y + p_rect.height < h_midpoint
        bot_quad = p_rect.y > h_midpoint

        if p_rect.x < v_midpoint and p_rect.x + p_rect.width < v_midpoint:
            if top_quad:
                index = 1
            elif bot_quad:
                index = 2
        elif p_rect.x > v_midpoint:
            if top_quad:
                index = 0
            elif bot_quad:
                index = 3

        return index

    def insert(self, p_rect: Rectangle):
        if self.nodes[0] is not None:
            index = self.get_index(p_rect)

            if index != -1:
                self.nodes[index].insert(p_rect)

                return

        self.objects.append(p_rect)

        if len(self.objects) > Quadtree.max_objects and self.level < Quadtree.max_levels:
            if self.nodes[0] is None:
                self.split()

            i = 0
            while i < len(self.objects):
                index = self.get_index(self.objects[i])
                if index != -1:
                    self.nodes[index].insert(self.objects.pop(i))
                else:
                    i += 1

    def retreive(self, return_obj: list, p_rect: Rectangle):
        index = self.get_index(p_rect)
        if index != -1 and self.nodes[0] is not None:
            self.nodes[index].retreive(return_obj, p_rect)

        return_obj.extend(self.objects)

        return return_obj

# quad = Quadtree(0, Rectangle(0, 0, 100, 100))


