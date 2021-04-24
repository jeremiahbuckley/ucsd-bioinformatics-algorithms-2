#! /bin/python3

import sys
import time
import random 

LEAD_FOLLOW_NODE_DELIMITER = " -> "
NODE_LIST_DELIMITER = ","


class NodeLeadFollowList:

    def __init__(self):
        self.leadingNodes = []
        self.followingNodes = []
'''
builds a dictionary like this:
[ACGT|GAGA] -> [AACGT|AGAGA] -> ["ACGTA|GAGAA", "ACGTA|GAGAC", .... "ACGTC|GAGAA", ..... ACGTG|GAGAA", etc...]
               [AACGT|CGAGA] -> [
               [AACGT|etc...
               [CACGT|AGAGA] -> ["ACGTA|GAGAA",.... "ACGTC|GAGAA"...., "ACGTC|GAGAC".... ]
               [CACGT|etc...
               [GACGT|AGAGA]etc..
               [TACGT|AGAGA]etc..

'''
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
                                    " all:" + NODE_LIST_DELIMITER.join(self.outgoing_nodes) + \
                                    " unused:" + NODE_LIST_DELIMITER.join(self.unused_outgoing_nodes)

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
                                    " all:" + NODE_LIST_DELIMITER.join(self.outgoing_nodes) + \
                                    " unused:" + NODE_LIST_DELIMITER.join(self.unused_outgoing_nodes) + \
                                    "  out:" + o
        return outbound

def node_lead_follow_pairs_to_list(node_lead_follow_lookup):
    out = []
    for kvp in node_lead_follow_lookup.items():
        #print(kvp[0])
        for edge_lead_node in kvp[1].leadingNodes:
            if len(kvp[1].followingNodes) > 0:
                out.append("{0} -> {1}".format(edge_lead_node, NODE_LIST_DELIMITER.join(kvp[1].followingNodes)))
    return out

def get_suffix_key_from_read_pair(read_pair, kmer_len):
    keys = read_pair.split("|")
    key = keys[0][1:kmer_len] + "|" + keys[1][1:kmer_len]
    return key

def get_prefix_key_from_read_pair(read_pair, kmer_len):
    keys = read_pair.split("|")
    key = keys[0][0:kmer_len-1] + "|" + keys[1][0:kmer_len-1]
    return key

def get_suffix_key_from_kmer(kmer, kmer_len):
    key = kmer[1:kmer_len]
    return key

def get_prefix_key_from_kmer(kmer, kmer_len):
    key = kmer[0:kmer_len-1]
    return key

def build_fragment_graph(fragments, fragments_are_read_pairs, kmer_len, read_distance):
    node_lead_follow_lookup = {}

    for i in range(len(fragments)):
        fragment = fragments[i]
        if fragments_are_read_pairs:
            fragment_suffix = get_suffix_key_from_read_pair(fragment, kmer_len)
            fragment_prefix = get_prefix_key_from_read_pair(fragment, kmer_len)
        else:
            fragment_suffix = get_suffix_key_from_kmer(fragment, kmer_len)
            fragment_prefix = get_prefix_key_from_kmer(fragment, kmer_len)

        if fragment_suffix not in node_lead_follow_lookup:
            node_lead_follow_lookup[fragment_suffix] = NodeLeadFollowList()
        if fragment not in node_lead_follow_lookup[fragment_suffix].leadingNodes:
            node_lead_follow_lookup[fragment_suffix].leadingNodes.append(fragment)

        if fragment_prefix not in node_lead_follow_lookup:
            node_lead_follow_lookup[fragment_prefix] = NodeLeadFollowList()
        if fragment not in node_lead_follow_lookup[fragment_prefix].followingNodes:
            node_lead_follow_lookup[fragment_prefix].followingNodes.append(fragment)

        #print(".." + fragment)
        #print(".." + fragment_suffix)
        #print(".." + fragment_prefix)
        #print_node_lead_follow_pairs(node_lead_follow_lookup)        
    return node_lead_follow_lookup



def init_data(nodes_list):
    nodes_data = {}

    for node_desc in nodes_list:
         #print(node_desc)
         end_id = node_desc.index(LEAD_FOLLOW_NODE_DELIMITER) #intentionally throw error if this isn't present
         nid = node_desc[0:end_id]
         outgoing_nodes = node_desc[end_id+4:len(node_desc)+1].rstrip()
         nodes_data[nid] = NodeInfo(nid, outgoing_nodes.split(NODE_LIST_DELIMITER))

    added_path = None
    for kvp in nodes_data.items():
        for out in kvp[1].outgoing_nodes:
            if out in nodes_data:
                nodes_data[out].incoming_nodes.append(kvp[0])
            else:
                added_path = (out, None)

    if added_path:
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
                                      NODE_LIST_DELIMITER.join(kvp[1].outgoing_nodes),
                                      NODE_LIST_DELIMITER.join(kvp[1].unused_outgoing_nodes)))


def rearrange_cycle(cycle_path, added_path):
    #print(added_path)
    new_path = cycle_path
    if added_path:
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

def find_eulerian_path(nodes_list, kmer_len, read_distance):
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

def reconstruct_dna_from_read_pair_list(read_pair_list, kmer_len, distance):
    read_pair = read_pair_list[0].split("|")
    strands = ["",""]
    for idx in [0,1]:
        strands[idx] = read_pair[idx]
    for i in range(1, len(read_pair_list)):
        read_pair = read_pair_list[i].split("|")
        for idx in [0,1]:
            strands[idx] += read_pair[idx][kmer_len-1:kmer_len]
    #print(strands[0])
    #print(strands[1])
    follow_strand_suffix = strands[1][len(strands[1])-kmer_len-distance:len(strands[1])]
    #print("      " + follow_strand_suffix)
    return strands[0] + follow_strand_suffix

def find_contigs_from_edges(edges_list, kmer_len):
    incoming_edge_count = {}
    for edge in edges_list:
        follow_nodes =edge[edge.index(LEAD_FOLLOW_NODE_DELIMITER)+4:len(edge)].split(NODE_LIST_DELIMITER)
        for follow_node in follow_nodes:
            if follow_node not in incoming_edge_count:
                incoming_edge_count[follow_node] = 0
            incoming_edge_count[follow_node] += 1

    print(incoming_edge_count)

    contigs_by_final_node = {}
    unclaimed_follow_nodes = []

    for edge in edges_list:
        lead_node = edge[0:edge.index(LEAD_FOLLOW_NODE_DELIMITER)]
        follow_node_list = edge[edge.index(LEAD_FOLLOW_NODE_DELIMITER) + 4: len(edge)]\
                          .split(NODE_LIST_DELIMITER)
        lead_node_starts_contig = True
        if lead_node in incoming_edge_count and incoming_edge_count[lead_node] == 1:
            lead_node_starts_contig = False

        print("_{0} -> {1}>".format(lead_node, "-".join(follow_node_list)))
        if lead_node_starts_contig:
            accounted_for_lead_node = False
            for follow_node in follow_node_list:
                if incoming_edge_count[follow_node] == 1:
                    contigs_by_final_node[follow_node] = [lead_node, final_node]
                    accounted_for_lead_node = True
                else:
                    if follow_node not in contigs_by_final_node:
                        contigs_by_final_node[follow_node] = [follow_node]
                        accounted_for_lead_node = True
                    else:
                        # this should happen occasisonaly, for nodes with > 1 incoming edge
                        pass
        else:
            if lead_node in contigs_by_final_node:
                for follow_node in follow_node_list:
                    if incoming_edge_count[follow_node] == 1:
                        nlist = contigs_by_final_node[lead_node]
                        nlist.append(follow_node_list[0])
                        del contigs_by_final_node[lead_node]
                        contigs_by_final_node[follow_node_list[0]] = nlist
                    else:
                        if follow_node not in contigs_by_final_node:
                            contigs_by_final_node[follow_node] = [follow_node]
                        else:
                            # this should happen occasisonaly, for nodes with > 1 incoming edge
                            pass
            else:
                raise ValueError("This is an unexpected place to end up. {0} -> {1}".format(lead_node, NODE_LIST_DELIMITER.join(follow_node_list)))
                     
    contigs = []
    for kvp in contigs_by_final_node:
        contig = kvp[1][0]
        if len(kvp[1]) > 1:
            for i in range(1,len(kvp[1])):
                contig += kvp[1][i][kmer_len-2,kmer_len]
        print(contig)
        contigs.append(contig)
    print("--")
    return contigs

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: one param, filename, format:<int_kmer_len> <int_pair_distance>\n<list_str_of_kmer_edges>.\nPositional switch -n will stop reading 1st line as kmer_len.")

    else:
        start = time.process_time()

        first_line_is_kmer_len = True
        if len(sys.argv) >= 3:
            if sys.argv[2] == "-n":
                first_line_is_kmer_len = False

        with open(sys.argv[1]) as f:
            edges = [line.rstrip() for line in f]

        distance = 1
        if not first_line_is_kmer_len:
            kmer_len = int(len(edges[0]) - 1 / 2)
        else:
            kmer_len = -1
            try:
                ints = [int(i) for i in edges[0].split(" ")]
                kmer_len = ints[0]
                distance = ints[1]
                edges.pop(0)
            except ValueError:
                kmer_len = int(len(edges[0])-1 / 2)

        node_lead_follow_structure = build_fragment_graph(edges, False, kmer_len, distance)
        lead_follow_list = node_lead_follow_pairs_to_list(node_lead_follow_structure)
        print(lead_follow_list)
        contigs_list = find_contigs_from_edges(lead_follow_list, kmer_len)
        for contig in contigs_list:
            print(contig)
        #read_pair_path = find_eulerian_path(lead_follow_list, kmer_len, distance)
        #print("->".join(read_pair_path))
        #full_dna = reconstruct_dna_from_read_pair_list(read_pair_path, kmer_len, distance)
        #print(full_dna)


        end = time.process_time()
        #print("Time: {0}".format(end-start))
