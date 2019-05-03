"""
IDDFS implemented in terms of a
recursive depth-limited DFS (called DLS).
"""
class vertex(self, id, attr):
	self.id = id
	self.attr = attr
	self.neighbours = set()
	def get_id(self):
		return self.id
	def get_attr(self):
		return self.attr
	def add_neighbour(self, neighbour):
		self.neighbours.add(neighbour)
	def get_neighbour(self):
		return self.neighbours
	def get_children(self, parent):
		return [neighbour for neighbour in self.neighbours if neighbour not in parent]


def IDDFS(root):
   for depth from 0 to ∞
       found, remaining ← DLS(root, depth, list())
       if found != None:
           return found
       elif not remaining:
           return None

def DLS(node, depth, visited):
   if depth = 0
       if node not in visited:
           return (node, True)
       else:
           return (None, True)    #(Not found, but may have children)

   elif depth > 0:
       any_remaining = False
       for child in node.get_children(): # to implement method of "children"
           found, remaining = DLS(child, depth−1, visited.append(node))
           if found != None
               return (found, True)
           if remaining:
               any_remaining = True    #(At least one node found at depth, let IDDFS deepen)
       return (None, any_remaining)
