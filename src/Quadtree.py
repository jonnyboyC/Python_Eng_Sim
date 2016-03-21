# Credit to Steven Lambert for the tutorial

class Rectangle:

    def __init__(self, x: float, y: float, width: float, height: float):
        self.x = float(x)
        self.y = float(y)
        self.width = float(width)
        self.height = float(height)


class Quadtree:

    max_objects = 20
    max_levels = 6

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

    def get_index(self, p_rect):
        index = -1
        bounds = p_rect.aabb
        v_midpoint = self.bounds.x + self.bounds.width/2
        h_midpoint = self.bounds.y + self.bounds.height/2

        top_quad = bounds.y < h_midpoint and bounds.y + bounds.height < h_midpoint
        bot_quad = bounds.y > h_midpoint

        if bounds.x < v_midpoint and bounds.x + bounds.width < v_midpoint:
            if top_quad:
                index = 1
            elif bot_quad:
                index = 2
        elif bounds.x > v_midpoint:
            if top_quad:
                index = 0
            elif bot_quad:
                index = 3

        return index

    def insert(self, p_rect):
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

    def retreive(self, return_obj: list, p_rect):
        index = self.get_index(p_rect)
        if index != -1 and self.nodes[0] is not None:
            self.nodes[index].retreive(return_obj, p_rect)

        return_obj.extend(self.objects)

        return return_obj

# quad = Quadtree(0, Rectangle(0, 0, 100, 100))


