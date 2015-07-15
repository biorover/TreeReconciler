#!/usr/bin/python
from ete2 import Tree
from math import log
import sys
import string




def count_ancestral_copynumber(tree,leaves_from_one_half,leaves_from_outgroups,pseudolist=[],return_nodes=False):
    count=0
    losscount_half1=0
    losscount_half2=0
    treeleaves=tree.get_leaves()
    node=treeleaves[0]
    passlist=list(set(treeleaves) & set(leaves_from_outgroups))
    leaves_other_half=list(set(treeleaves)-set(leaves_from_one_half)-set(leaves_from_outgroups))
    pseudoleaves=[]
    nodelist=[]
    for i in treeleaves:
        if i.name in pseudolist:
            pseudoleaves.append(i)
    non_pseudo_leaves_half1=list(set(leaves_from_one_half)-set(pseudoleaves))
    non_pseudo_leaves_half2=list(set(leaves_other_half)-set(pseudoleaves))
    while len(passlist) < len(treeleaves):
        sislist=node.get_sisters()
        sisleaves=[]
        for i in sislist:
            sisleaves=sisleaves+i.get_leaves()
        nodeleaves=node.get_leaves()
        #all children in pass list
        if set(nodeleaves) <= set(passlist):
            try:
                node=list(set(treeleaves) - set(passlist))[0]
            except:
                pass
        #All children in B1
        elif set(nodeleaves) <= set(leaves_from_one_half):
            #All sisters in B1
            if set(sisleaves) <= set(leaves_from_one_half):
                try: node=node.get_ancestors()[0]
                except IndexError:
                    count=count+1
                    break
            #Some sisters in B2 or outgroup: determine which sisters not statistically different from self
            else:
                count=count+1
                nodelist.append(Tree())
                nodelist[-1].add_child(node.copy())
                passlist=passlist+nodeleaves
                loss_half2=True
                for i in sislist:
                    ileaves=i.get_leaves()
                    if set(ileaves) <= set(leaves_from_one_half):
                        passlist=passlist+ileaves
                        nodelist[-1].add_child(i.copy())
                    elif set(ileaves) <= set(leaves_other_half):
                        passlist=passlist+ileaves
                        nodelist[-1].add_child(i.copy())
                        loss_half2=False
                if loss_half2:
                    losscount_half2=losscount_half2+1
                if set(nodeleaves)-set(pseudoleaves)==set():
                    losscount_half1=losscount_half1+1
                try:
                    node=list(set(treeleaves) - set(passlist))[0]
                except:
                    pass
                if len(nodelist[-1].children) == 1:
                    nodelist[-1]=nodelist[-1].children[0]
        #All children in B2
        elif set(nodeleaves) <= set(leaves_other_half):
            #All sisters in B2
            if set(sisleaves) <= set(leaves_other_half):
                try: node=node.get_ancestors()[0]
                except IndexError:
                    count=count+1
                    break
            #Some sisters in B1 or outgroup: determine which sisters not statistically different from self
            else:
                count=count+1
                nodelist.append(Tree())
                nodelist[-1].add_child(node.copy())
                passlist=passlist+nodeleaves
                loss_half1=True
                for i in sislist:
                    ileaves=i.get_leaves()
                    if set(ileaves) <= set(leaves_other_half):
                        passlist=passlist+ileaves
                        nodelist[-1].add_child(i.copy())
                    elif set(ileaves) <= set(leaves_from_one_half):
                        passlist=passlist+ileaves
                        nodelist[-1].add_child(i.copy())
                        loss_half1=False
                if loss_half1:
                    losscount_half1=losscount_half1+1
                if set(nodeleaves)-set(pseudoleaves)==set():
                    losscount_half2=losscount_half2+1
                try:
                    node=list(set(treeleaves) - set(passlist))[0]
                except:
                    pass
                if len(nodelist[-1].children) == 1:
                    nodelist[-1]=nodelist[-1].children[0]
        #Children from B1 and B2 or either + outgroup
        else:
            count=count+1
            focal_leaves=set(nodeleaves)-set(leaves_from_outgroups)
            prunelist=[]
            for i in focal_leaves:
                prunelist.append(i.name)
            nodelist.append(node.copy().prune(prunelist))
            if set(nodeleaves) <= set(leaves_from_outgroups) | set(non_pseudo_leaves_half1):
                losscount_half2=losscount_half2+1
            elif set(nodeleaves) <= set(leaves_from_outgroups) | set(non_pseudo_leaves_half2):
                losscount_half1=losscount_half1+1
            passlist=passlist+nodeleaves
    if not return_nodes:
        return (count,losscount_half1,losscount_half2)
    else:
        return nodelist


def birthdeath_rates(species_tree,gene_tree,gene_to_species_map,pseudolist=[],clade_focalnode=0):
    #build species genelists
    speciesdic={}
    nonpseudodic={}
    for i in species_tree.get_leaves():
        speciesdic[i.name]=[]
        nonpseudodic[i.name]=[]
    for i in gene_to_species_map:
        try:
            speciesdic[i[0]].append(gene_tree.get_leaves_by_name(i[1])[0])
            if not i[1] in pseudolist:
                nonpseudodic[i[0]].append(i[1])
        except:
            pass
    #build species split database
    split_names=[]
    split_distances=[]
    split_copynumbers=[]
    split_birthrates=[]
    split_deathrates=[]
    split_BDrates=[]
    split_births=[]
    split_deaths=[]
    doubles=[]
    halves=[]
    nodelist=[species_tree]+species_tree.get_descendants()
    nodedic={}
    for node in nodelist:
        if not node.is_leaf():
            b1=node.children[0]
            b2=node.children[1]
            nodedic[node]=[len(split_names),b1,b2]
            og=list(set(species_tree.get_leaves()).difference(set(b1.get_leaves()) | set (b2.get_leaves())))
            split_names.append((b1.get_leaf_names(),b2.get_leaf_names()))
            split_distances.append(node.dist)
            halfleaves=[]
            ogleaves=[]
            for taxon in b1.get_leaf_names():
                halfleaves=halfleaves+speciesdic[taxon]
            for taxon in og:
                ogleaves=ogleaves+speciesdic[taxon.name]
            split_copynumbers.append(count_ancestral_copynumber(gene_tree,halfleaves,ogleaves,pseudolist=pseudolist))
            if len(split_names) > 1:
                parent_data=nodedic[node.up]
                parent_index=parent_data[0]
                if node == parent_data[1]:
                    bindi=1
                elif node == parent_data[2]:
                    bindi=2
                else:
                    print "didn't work"
        else:
            split_names.append([node.name])
            split_distances.append(node.dist)
            split_copynumbers.append((len(nonpseudodic[node.name]),0,0))
            parent_data=nodedic[node.up]
            parent_index=parent_data[0]
            if node == parent_data[1]:
                bindi=1
            elif node == parent_data[2]:
                bindi=2
        if len(split_names) > 1:
            anCpN=split_copynumbers[parent_index][0]-split_copynumbers[parent_index][bindi]
            branchCpN=split_copynumbers[-1][0]
            deltCpN=branchCpN-anCpN
            split_births.append(deltCpN)
            split_deaths.append(split_copynumbers[parent_index][0]-anCpN)
            if anCpN > 0 and branchCpN > 0:
                split_birthrates.append(log((branchCpN*1.0/anCpN),2)/split_distances[-1])
                doubles.append(log((branchCpN*1.0/anCpN),2))
            else:
                split_birthrates.append(0.0)
                doubles.append(0.0)
            if anCpN > 0:
                split_deathrates.append(abs(log((anCpN*1.0/split_copynumbers[parent_index][0]),2)/split_distances[-1]))
                halves.append(abs(log((anCpN*1.0/split_copynumbers[parent_index][0]),2)))
            elif split_copynumbers[parent_index][0] > 0:
                split_deathrates.append((1.0-log(1.0/split_copynumbers[parent_index][0],2))/split_distances[-1])
                halves.append(1.0-log(1.0/split_copynumbers[parent_index][0],2))
            else:
                split_deathrates.append(0.0)
                halves.append(0.0)
            if branchCpN > 0:
                split_BDrates.append(log((branchCpN*1.0/split_copynumbers[parent_index][0]),2)/split_distances[-1])
            elif split_copynumbers[parent_index][0] > 0:
                split_BDrates.append(-1*(1.0-log(1.0/split_copynumbers[parent_index][0],2))/split_distances[-1])
            else:
                split_BDrates.append(0.0)
    totalbirths=sum(split_births[clade_focalnode:])
    totaldeaths=sum(split_deaths[clade_focalnode:])
    totaldoubles=sum(doubles[clade_focalnode:])
    totalhalves=sum(halves[clade_focalnode:])
    totaldistance=sum(split_distances[clade_focalnode+1:])
    totalanCpN=split_copynumbers[clade_focalnode][0]
    clade_birthrate=totaldoubles/totaldistance
    clade_deathrate=totalhalves/totaldistance
    return [split_names,split_distances,split_copynumbers,split_birthrates,split_deathrates,split_BDrates,doubles,halves,clade_birthrate,clade_deathrate]








def write2nexus(ete2_tree,annotation_names,support=False):
    nexstring="#NEXUS\nbegin taxa;\n\tdimensions ntax="
    taxnames=ete2_tree.get_leaf_names()
    nexstring=nexstring+str(len(taxnames))+";\n\ttaxlabels\n\t"+"\n\t".join(taxnames)+"\n;\nend;\n\nbegin trees;\n\ttree tree_1 = [&R] "
    modtree=ete2_tree.copy()
    for node in modtree.get_descendants():
        anlist=[]
        for annotation in annotation_names:
            try:
                anlist.append(annotation+"equals_sign"+str(getattr(node,annotation)))
            except AttributeError:
                pass
        if node.name == "NoName":
            if anlist != []:
                node.name=str(anlist).replace("'","").replace('[','leftbracket').replace(']','rightbracket').replace(',','comma_mark')
            else:
                node.name=""
        else:
            if anlist != []:
                node.name=node.name+str(anlist).replace("'","").replace('[','leftbracket').replace(']','rightbracket').replace(',','comma_mark')
            else:
                pass
    nexstring=nexstring+modtree.write(format=1).replace('leftbracket','[&').replace('rightbracket',']').replace('equals_sign','=').replace('comma_mark',',')+"\nend;"
    return nexstring



def namecode2map(namecodelistortree,genetree):
    ormap=[]
    try:
        x=namecodelistortree.get_leaf_names()
        ncl=[]
        for i in x:
            ncl.append([i,i])
    except:
        ncl=namecodelistortree
    for i in genetree.get_leaves():
        for k in ncl:
            if k[0] in i.name:
                ormap.append([k[-1],i.name])
    return ormap
