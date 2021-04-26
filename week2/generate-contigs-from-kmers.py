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
                out.append("{0}{1}{2}".format(edge_lead_node, LEAD_FOLLOW_NODE_DELIMITER, NODE_LIST_DELIMITER.join(kvp[1].followingNodes)))
            else:
                out.append("{0}{1}{2}".format(edge_lead_node, LEAD_FOLLOW_NODE_DELIMITER,"None"))
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

def build_fragment_graph_old(fragments, fragments_are_read_pairs, kmer_len, read_distance):
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
        #print("...." + fragment_suffix)
        #print("......" + fragment_prefix)
    
    #for kvp in node_lead_follow_lookup.items():
        #print("{0} {1} ->{2}".format(kvp[0], kvp[1].leadingNodes, kvp[1].followingNodes))
    return node_lead_follow_lookup

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
        node_lead_follow_lookup[fragment_suffix].leadingNodes.append(fragment)

        if fragment_prefix not in node_lead_follow_lookup:
            node_lead_follow_lookup[fragment_prefix] = NodeLeadFollowList()
        node_lead_follow_lookup[fragment_prefix].followingNodes.append(fragment)

        #print(".." + fragment)
        #print("...." + fragment_suffix)
        #print("......" + fragment_prefix)
    
    #for kvp in node_lead_follow_lookup.items():
    #    print("{0} {1} ->{2}".format(kvp[0], kvp[1].leadingNodes, kvp[1].followingNodes))
    return node_lead_follow_lookup

def init_data(nodes_list):
    nodes_data = {}

    for node_desc in nodes_list:
         #print(node_desc)
         end_id = node_desc.index(LEAD_FOLLOW_NODE_DELIMITER) #intentionally throw error if this isn't present
         nid = node_desc[0:end_id]
         outgoing_nodes = node_desc[end_id+4:len(node_desc)+1].rstrip()
         if outgoing_nodes == "None":
             nodes_data[nid] = NodeInfo(nid, [])
         else:
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

def find_contigs_from_lead_follow_lookup(lead_follow_lookup, kmer_len):
    contigs_nodes = []
    possible_isolated_nodes = set()
    for kvp in lead_follow_lookup.items():
        node = kvp[1]
        if len(node.leadingNodes) <= 1 and len(node.followingNodes) <= 1:
            possible_isolated_nodes.add(kvp[0])

    for kvp in lead_follow_lookup.items():
        node = kvp[1]
        if len(node.leadingNodes) != 1 or len(node.followingNodes) != 1:
            if len(node.followingNodes) > 0:
                for next_node in node.followingNodes:
                    next_node_suffix_key = get_suffix_key_from_kmer(next_node, kmer_len)
                    non_branching_path = [kvp[0], next_node_suffix_key]
                    if kvp[0] in possible_isolated_nodes:
                        possible_isolated_nodes.remove(kvp[0])
                    if next_node_suffix_key in possible_isolated_nodes:
                        possible_isolated_nodes.remove(next_node_suffix_key)
                    nn = lead_follow_lookup[next_node_suffix_key]
                    while len(nn.leadingNodes) == 1 and len(nn.followingNodes) == 1:
                        nn_follow = nn.followingNodes[0]
                        next_node_suffix_key = get_suffix_key_from_kmer(nn_follow, kmer_len)
                        non_branching_path.append(next_node_suffix_key)
                        if next_node_suffix_key in possible_isolated_nodes:
                            possible_isolated_nodes.remove(next_node_suffix_key)
                        nn = lead_follow_lookup[next_node_suffix_key]
                    contigs_nodes.append(non_branching_path)

    #print(contigs_nodes)
    #print(possible_isolated_nodes)
    while len(possible_isolated_nodes) > 0:
        nodes_to_remove = []

        key = possible_isolated_nodes.pop()
        nodes_to_remove.append(key)
        next_node = lead_follow_lookup[key].followingNodes[0]
        next_node_suffix_key = get_suffix_key_from_kmer(next_node, kmer_len)
        if next_node_suffix_key in possible_isolated_nodes:
            nodes_to_remove.append(next_node_suffix_key)
        non_branching_path = [node, next_node_suffix_key]

        while next_node_suffix_key in possible_isolated_nodes:
            nn = lead_follow_lookup[next_node_suffix_key]
            nn_follow = nn.followingNodes[0]
            next_node_suffix_key = get_suffix_key_from_kmer(nn_follow, kmer_len)
            non_branching_path.append(next_node_suffix_key)
            if next_node_suffix_key in possible_isolated_nodes:
                nodes_to_remove.append(next_node_suffix_key)
            nn = lead_follow_lookup[next_node_suffix_key]
        
        contigs_nodes.append(non_branching_path)    

        for n in nodes_to_remove:
            possible_isolated_nodes.remove(n)

    #print(contigs_nodes)
    #print("--")
    contigs = []
    for kmers in contigs_nodes:
        #print(type(kmers))
        contig = kmers[0]
        if len(kmers) > 1:
            for i in range(1,len(kmers)):
                contig += kmers[i][kmer_len-2:kmer_len]
        contigs.append(contig)
        #print(contig)
    return contigs
        



def find_contigs_from_edges(edges_list, kmer_len):
    incoming_edge_count = {}
    for edge in edges_list:
        follow_nodes_str =edge[edge.index(LEAD_FOLLOW_NODE_DELIMITER)+4:len(edge)]
        follow_nodes = []
        if follow_nodes_str != "None":
            follow_nodes = follow_nodes_str.split(NODE_LIST_DELIMITER)
        for follow_node in follow_nodes:
            if follow_node not in incoming_edge_count:
                incoming_edge_count[follow_node] = 0
            incoming_edge_count[follow_node] += 1

    #print(incoming_edge_count)

    contigs_by_final_node = {}
    unclaimed_follow_nodes = []

    print("--")
    print(edges_list)
    contigs_nodes = []

    watch_node = "CGTCCAAGGCGGGTAGAATATCAAGCGGATGCTTACTAATCATCGAACCGTGGCGCAGCCGGTAGACG"

    possible_follow_nodes_that_need_to_be_added = []
    for edge in edges_list:
        lead_node = edge[0:edge.index(LEAD_FOLLOW_NODE_DELIMITER)]
        follow_node_str = edge[edge.index(LEAD_FOLLOW_NODE_DELIMITER) + 4: len(edge)]
        follow_node_list = []
        if follow_node_str != "None":
            follow_node_list = edge[edge.index(LEAD_FOLLOW_NODE_DELIMITER) + 4: len(edge)]\
                                .split(NODE_LIST_DELIMITER)

        if lead_node == watch_node or watch_node in follow_node_list:
            print("^" + edge)
            for contig in contigs_nodes:
                print("  {0}".format(",".join(contig)))

        if lead_node in incoming_edge_count and incoming_edge_count[lead_node] > 1:
           if len(follow_node_list) == 0:
               contigs_nodes.append([lead_node])
           else:
               added_one = False
               for follow_node in follow_node_list:
                   if incoming_edge_count[follow_node] > 1:
                       possible_follow_nodes_that_need_to_be_added.append(follow_node)
                       if not added_one:
                           contigs_nodes.append([lead_node])
                           added_one = True
                   else:
                       contigs_nodes.append([lead_node, follow_node])
                       added_one = True

        if lead_node == watch_node or watch_node in follow_node_list:
            for contig in contigs_nodes:
                print("      {0}".format(",".join(contig)))

    for node in possible_follow_nodes_that_need_to_be_added:
        not_found = True
        for contig in contigs_nodes:
            if contig[0] == node:
                not_found = False
        if not_found:
            contigs_nodes.append([node])

    for edge in edges_list:
        lead_node = edge[0:edge.index(LEAD_FOLLOW_NODE_DELIMITER)]
        follow_node_str = edge[edge.index(LEAD_FOLLOW_NODE_DELIMITER) + 4: len(edge)]
        follow_node_list = []
        if follow_node_str != "None":
            follow_node_list = edge[edge.index(LEAD_FOLLOW_NODE_DELIMITER) + 4: len(edge)]\
                                .split(NODE_LIST_DELIMITER)
        if lead_node == watch_node or watch_node in follow_node_list:
            print("^" + edge)
            for contig in contigs_nodes:
                print("  {0}".format(",".join(contig)))

        if lead_node in incoming_edge_count and incoming_edge_count[lead_node] > 1:
            pass
        else:
            contigs_w_lead_node_start = []
            contigs_w_lead_node_end = []
            for contig_idx in range(len(contigs_nodes)):
                if contigs_nodes[contig_idx][0] == lead_node:
                    contigs_w_lead_node_start.append(contig_idx)
                if contigs_nodes[contig_idx][len(contigs_nodes[contig_idx])-1] == lead_node:
                    contigs_w_lead_node_end.append(contig_idx)

            if len(follow_node_list) == 0:
                if len(contigs_w_lead_node_start) == 0 and len(contigs_w_lead_node_end) == 0:
                    contigs_nodes.append([lead_node])
                    contigs_w_lead_node_start.append(len(contigs_nodes))
                    contigs_w_lead_node_end.append(len(contigs_nodes))
            else:
                for follow_node in follow_node_list:
                    contigs_w_follow_node_start = []
                    contigs_w_follow_node_end = []
                    for contig_idx in range(len(contigs_nodes)):
                        if contigs_nodes[contig_idx][0] == follow_node:
                            contigs_w_follow_node_start.append(contig_idx)
                        if contigs_nodes[contig_idx][len(contigs_nodes[contig_idx])-1] == follow_node:
                            contigs_w_follow_node_end.append(contig_idx)
                    follow_node_has_multiple_inroutes = incoming_edge_count[follow_node] > 1

                    for i in contigs_w_follow_node_end:
                        if i not in contigs_w_follow_node_start:
                            if follow_node_has_multiple_inroutes:
                                raise RuntimeError("Found contig size > 1 with multi-input-at-edge. kmer = {0}, contig = {1}".format(follow_node, contigs_nodes[i]))
                            else:
                                raise RuntimeError("Found second inroute edge to {0} from {1}. Already existing is: {2}".format(follow_node, lead_node, contigs_nodes[i]))

                    if follow_node_has_multiple_inroutes:
                        if len(contigs_w_follow_node_start) > 0:
                           if len(contigs_w_lead_node_start) == 0 and len(contigs_w_lead_node_end) == 0:
                               contigs_nodes.append([lead_node])
                               contigs_w_lead_node_start.append(len(contigs_nodes)-1)
                               contigs_w_lead_node_end.append(len(contigs_nodes)-1)
                        else:
                            contigs_nodes.append([follow_node])
                            contigs_w_follow_node_start.append(len(contigs_nodes)-1)
                            contigs_w_follow_node_end.append(len(contigs_nodes)-1)
                            if len(contigs_w_lead_node_start) == 0 and len(contigs_w_lead_node_end) == 0:
                               contigs_nodes.append([lead_node])
                               contigs_w_lead_node_start.append(len(contigs_nodes)-1)
                               contigs_w_lead_node_end.append(len(contigs_nodes)-1)
                    else:
                        if len(contigs_w_follow_node_start) == 0:
                            if len(contigs_w_lead_node_end) == 0:
                                contigs_nodes.append([lead_node, follow_node])
                                contigs_w_lead_node_start.append(len(lead_node))
                                contigs_w_follow_node_end.append(len(follow_node))
                            else:
                                idxs_to_remove = []
                                for idx in contigs_w_lead_node_end:
                                    contigs_nodes[idx].append(follow_node)
                                    idxs_to_remove.append(idx)
                                    contigs_w_follow_node_end.append(idx)
                                for idx in idxs_to_remove:
                                    contigs_w_lead_node_end.remove(idx)
                        else:
                            # for now, do "insert at beginning loop" separately from "join two contigs" loop, just to confirm logic is good
                            for i in contigs_w_follow_node_start:
                                idxs_to_remove = []
                                if i not in contigs_w_lead_node_end:
                                    contigs_nodes[i].insert(0, lead_node)
                                    idxs_to_remove.append(i)
                                    contigs_w_lead_node_start.append(i)
                                for idx in idxs_to_remove:
                                    contigs_w_follow_node_start.remove(idx)
                            for i in contigs_w_follow_node_start:
                                idxs_to_remove = []
                                if i in contigs_w_lead_node_end:
                                   ln_idx = contigs_w_lead_node_end(i)
                                   if (len(contigs_nodes[i]) > 1):
                                       contigs_nodes[ln_idx].extend(contigs_nodes[i][1:len(contigs_nodes[i])])
                                   idxs_to_remove.append(i)
                                   contigs_w_lead_node_end.pop(ln_idx)
                                for idx in idxs_to_remove:
                                   contigs_w_follow_node_start.remove(idx)

        if lead_node == watch_node or watch_node in follow_node_list:
            for contig in contigs_nodes:
                print("    {0}".format(",".join(contig)))

    #print("--")
    contigs = []
    for kmers in contigs_nodes:
        #print(type(kmers))
        contig = kmers[0]
        if len(kmers) > 1:
            for i in range(1,len(kmers)):
                contig += kmers[i][kmer_len-1:kmer_len]
        contigs.append(contig)
        #print(contig)
    return contigs
'''
    while len(unclassified_edges_list) > 0:
        still_unclassified_edges = []
        for edge in unclassified_edges_list:
            lead_node = edge[0:edge.index(LEAD_FOLLOW_NODE_DELIMITER)]
            follow_node_str = edge[edge.index(LEAD_FOLLOW_NODE_DELIMITER) + 4: len(edge)]
            follow_node_list = []
            if follow_node_str != "None":
                follow_node_list = edge[edge.index(LEAD_FOLLOW_NODE_DELIMITER) + 4: len(edge)]\
                                    .split(NODE_LIST_DELIMITER)
            lead_node_starts_contig = True
            if lead_node in incoming_edge_count and incoming_edge_count[lead_node] == 1:
                lead_node_starts_contig = False

            print("_{0} -> {1}".format(lead_node, ",".join(follow_node_list)))
            if lead_node_starts_contig:
                accounted_for_lead_node = False
                for follow_node in follow_node_list:
                    if incoming_edge_count[follow_node] == 1:
                        contigs_by_final_node[follow_node] = [lead_node, follow_node]
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
                    print(follow_node_list)
                    for follow_node in follow_node_list:
                        print(follow_node)
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
                    still_unclassified_edges.append(edge)

        if len(still_unclassified_edges) > 0:
            print("\n**  ".join(still_unclassified_edges))
        unclassified_edges_list = still_unclassified_edges
'''
                     
    #print(contigs_by_final_node)
    #contigs = []
    #for kvp in contigs_by_final_node.items():
    #    kmers = kvp[1]
    #    print(type(kmers))
    #    contig = kmers[0]
    #    if len(kmers) > 1:
    #        for i in range(1,len(kmers)):
    #            contig += kmers[i][kmer_len-1:kmer_len]
    #    print(contig)
    #    contigs.append(contig)
    #print("--")
    #return contigs

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
            kmer_len = int(len(edges[0]))
        else:
            kmer_len = -1
            try:
                ints = [int(i) for i in edges[0].split(" ")]
                kmer_len = ints[0]
                distance = ints[1]
                edges.pop(0)
            except ValueError:
                kmer_len = int(len(edges[0]))

        node_lead_follow_structure = build_fragment_graph(edges, False, kmer_len, distance)
        contigs_list = find_contigs_from_lead_follow_lookup(node_lead_follow_structure, kmer_len)
        #lead_follow_list = node_lead_follow_pairs_to_list(node_lead_follow_structure)
        #print(lead_follow_list)
        #print("x")
        #contigs_list = find_contigs_from_edges(lead_follow_list, kmer_len)
        for contig in contigs_list:
            print(contig)
        #read_pair_path = find_eulerian_path(lead_follow_list, kmer_len, distance)
        #print("->".join(read_pair_path))
        #full_dna = reconstruct_dna_from_read_pair_list(read_pair_path, kmer_len, distance)
        #print(full_dna)


        end = time.process_time()
        #print("Time: {0}".format(end-start))
