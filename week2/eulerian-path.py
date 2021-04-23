#! /bin/python3

import sys
import time
import random

class NodeInfo:
    def __init__(self, nid, outgoing_nodes):
        self.node_id = nid
        self.outgoing_nodes = outgoing_nodes[:]
        self.unused_outgoing_nodes = outgoing_nodes[:]
        self.incoming_nodes = []
        self.last_cycle_debug_str = ""

    def has_open_edges(self):
        return len(self.unused_outgoing_nodes) > 0

    def is_balanced(self):
        return len(self.outgoing_nodes) == len(self.incoming_nodes)

    def next_trip_out(self):
        self.last_cycle_debug_str = "   " + \
                                    " id:" + self.node_id + \
                                    " all:" + ",".join(self.outgoing_nodes) + \
                                    " unused:" + ",".join(self.unused_outgoing_nodes)

        outbound = None
        o_idx = 0
        if len(self.unused_outgoing_nodes) > 1:
            o_idx = random.randrange(len(self.unused_outgoing_nodes))
        outbound = self.unused_outgoing_nodes.pop(o_idx)

        o = outbound
        if not outbound:
            o = "None"
        self.last_cycle_debug_str+= "..." + \
                                    " id:" + self.node_id + \
                                    " all:" + ",".join(self.outgoing_nodes) + \
                                    " unused:" + ",".join(self.unused_outgoing_nodes) + \
                                    "  out:" + o
        return outbound

def init_data(nodes_list):
    nodes_data = {}

    for node_desc in nodes_list:
         end_id = node_desc.index(" -> ") #intentionally throw error if this isn't present
         nid = node_desc[0:end_id]
         outgoing_nodes = node_desc[end_id+4:len(node_desc)+1].rstrip()
         nodes_data[nid] = NodeInfo(nid, outgoing_nodes.split(","))

    added_path = ()
    for kvp in nodes_data.items():
        for out in kvp[1].outgoing_nodes:
            if out in nodes_data:
                nodes_data[out].incoming_nodes.append(kvp[0])
            else:
                added_path = (out, None)

    for kvp in nodes_data.items():
        if len(kvp[1].outgoing_nodes) != len(kvp[1].incoming_nodes):
            added_path = (added_path[0], kvp[0])
            break;

    nodes_data[added_path[0]] = NodeInfo(added_path[0], [added_path[1]])


    #print_state(nodes_data, 0)

    return nodes_data, added_path

def print_state(nodes_data, indent):

    for kvp in nodes_data.items():
        print("{0}{1}-{2}-{3}".format(" " * indent,
                                      kvp[0],
                                      ",".join(kvp[1].outgoing_nodes),
                                      ",".join(kvp[1].unused_outgoing_nodes)))


def rearrange_cycle(cycle_path, added_path):
    #print(added_path)
    out_idx = cycle_path.index(added_path[0])
    new_path = cycle_path[out_idx + 1:len(cycle_path)] + \
               cycle_path[1:out_idx+1]
    return new_path

def restart_cycle(nodes_data, path):
    can_restart = False
    restart_node_key = None

    restart_nodes = []

    for node_key in path:
        if nodes_data[node_key].has_open_edges():
            restart_nodes.append(node_key)
            can_restart = True

    #for nk in restart_nodes:
    #    print(nodes_data[nk].last_cycle_debug_str)

    if len(restart_nodes) > 0:
        restart_point = 0
        if len(restart_nodes) > 1:
            restart_point = random.randrange(len(restart_nodes))
        restart_node_key = restart_nodes[restart_point]

    new_path = path
    if can_restart:
        # intetionally cutting out last char, since it's also the first char
        # the split-and-rejoin would mean a duplicate at the join point
        new_path = path[path.index(restart_node_key):len(path)-1]+ \
                   path[0:path.index(restart_node_key)]

    #print("restarting {0} {1}".format(can_restart, restart_node_key))
    return can_restart, restart_node_key, new_path

def find_eulerian_path(nodes_list):
    nodes_data, added_path = init_data(nodes_list)

    cycle_path = []
    start_idx = random.randrange(len(nodes_data))
    old_node_key = None
    for kvp in nodes_data.items():
        if start_idx == 0:
            old_node_key = kvp[0]
        start_idx -= 1

    graph_has_open_edges = nodes_data[old_node_key].has_open_edges()
    cycle_path.append(old_node_key)
    while graph_has_open_edges:
        #print(cycle_path)
        #print_state(nodes_data, len(cycle_path))
        if nodes_data[old_node_key].has_open_edges():
            new_node_key = nodes_data[old_node_key].next_trip_out()
            cycle_path.append(new_node_key)
        else:
            #print("search needs restart-" + old_node_key + "\n" + \
            #      nodes_data[old_node_key].last_cycle_debug_str)
            graph_has_open_edges, new_node_key, cycle_path = restart_cycle(nodes_data, cycle_path)
            if graph_has_open_edges:
                cycle_path.append(new_node_key)
 
        old_node_key = new_node_key

    cycle_path = rearrange_cycle(cycle_path, added_path)

    return cycle_path


if __name__ == "__main__":
    start = time.process_time()
   
    random.seed(0)

    if len(sys.argv) < 2:
        print("Usage: expect file input of format: <string_list_of_nodes>.\n" \
              "Node format is [parent] -> [child1],[child2],...,[childN]")

    with open(sys.argv[1]) as fl:
        nodes = fl.readlines()

    #print(nodes)

    path = find_eulerian_path(nodes)
    print("->".join(path))

    end = time.process_time()
    #print("Time: {0}".format(end-start))
