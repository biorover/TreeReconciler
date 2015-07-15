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



if __name__ == "__main__":
    gT=Tree(open('RAxML_Threshold-70-ConsensusTree.rbs70.consensus.tre').read())
    litdated_sptree=Tree(open('Ronquist_Moreau_dates.tre').read())
    orcosubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr1')[0],gT.get_leaves_by_name('_seed_PbarOr1')[0]])
    gT.set_outgroup(orcosubfam)
    csubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_AmelOr116')[0],gT.get_leaves_by_name('_seed_PbarOr197')[0]])
    dsubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr122')[0],gT.get_leaves_by_name('_seed_PbarOr229')[0]])
    rsubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_HsalOr179')[0],gT.get_leaves_by_name('_seed_PbarOr148')[0]])
    ssubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr291')[0],gT.get_leaves_by_name('_seed_PbarOr149')[0]])
    qsubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_AmelOr160')[0],gT.get_leaves_by_name('_seed_PbarOr146')[0]])
    ksubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr2')[0],gT.get_leaves_by_name('_seed_PbarOr3')[0]])
    lsubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr5')[0],gT.get_leaves_by_name('_seed_PbarOr39')[0]])
    tsubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr20')[0],gT.get_leaves_by_name('_seed_PbarOr72NTE')[0]])
    usubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr38')[0],gT.get_leaves_by_name('_seed_PbarOr92NTE')[0]])
    msubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_HsalOr60')[0],gT.get_leaves_by_name('_seed_PbarOr53')[0]])
    nsubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr61')[0],gT.get_leaves_by_name('_seed_PbarOr65')[0]])
    psubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_AmelOr63')[0],gT.get_leaves_by_name('_seed_PbarOr63')[0]])
    osubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr41')[0],gT.get_leaves_by_name('_seed_CbirOr50')[0]])
    vsubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr62')[0],gT.get_leaves_by_name('_seed_PbarOr104')[0]])
    asubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_LhumOr59')[0],gT.get_leaves_by_name('_seed_PbarOr175')[0]])
    jsubfam=gT.get_common_ancestor([gT.get_leaves_by_name('NvitOr74PSE')[0],gT.get_leaves_by_name('_seed_PbarOr154')[0]])
    isubfam=gT.get_common_ancestor([gT.get_leaves_by_name('NvitOr297PSE')[0],gT.get_leaves_by_name('_seed_PbarOr145')[0]])
    bsubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_CbirOr159')[0],gT.get_leaves_by_name('_seed_PbarOr151')[0]])
    gsubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr45')[0],gT.get_leaves_by_name('_seed_PbarOr181')[0]])
    hsubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr44')[0],gT.get_leaves_by_name('_seed_PbarOr178')[0]])
    nesubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr300')[0],gT.get_leaves_by_name('_seed_PbarOr268')[0]])
    fsubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr261')[0],gT.get_leaves_by_name('_seed_PbarOr162')[0]])
    esubfam=gT.get_common_ancestor([gT.get_leaves_by_name('_seed_NvitOr98')[0],gT.get_leaves_by_name('_seed_PbarOr201')[0]])
    amleaves=[]
    nvleaves=[]
    lhleaves=[]
    cfleaves=[]
    pbleaves=[]
    cbleaves=[]
    hsleaves=[]
    ormap=[]
    for i in gT.get_leaves():
     if 'Am' in i.name:
         amleaves.append(i)
         ormap.append(['Amel',i.name])
     elif 'Nv' in i.name:
         nvleaves.append(i)
         ormap.append(['Nvit',i.name])
     elif 'Lh' in i.name:
         lhleaves.append(i)
         ormap.append(['Lhum',i.name])
     elif 'Cf' in i.name:
         cfleaves.append(i)
         ormap.append(['Cflo',i.name])
     elif 'Pb' in i.name:
         pbleaves.append(i)
         ormap.append(['Pbar',i.name])
     elif 'Cb' in i.name:
         cbleaves.append(i)
         ormap.append(['Cbir',i.name])
     elif 'Hs' in i.name:
         hsleaves.append(i)
         ormap.append(['Hsal',i.name])
    sftreelist=[]
    subfamnames=[]
    for i in string.lowercase[:22]:
        sftreelist.append(eval(i+'subfam'))
        subfamnames.append(i)
    sftreelist.append(nesubfam)
    subfamnames.append('ne')
#birtharray=[]
#namesarray=[]
    ortab=open('/Volumes/antqueen/genomics/experiments/analyses/SKM20141001_ORevolution/AllOrsTab.txt').readlines()
    pseudolist=[]
    for i in ortab:
        if i.split('\t')[8]=="yes":
            pseudolist.append(i.split('\t')[0])
    allarray=[]
    for i in sftreelist:
        x=Get_ancestral_copy_numbers.species_branch_birthdeath_rates(litdated_sptree,i,ormap,pseudolist=pseudolist)
#  namesarray.append(x[0])
#  birtharray.append(x[3])
        allarray.append(x)
    antree=litdated_sptree.copy()
    for node in antree.get_descendants():
        nodeindi=antree.get_descendants().index(node)
        #print node
        #print nodeindi
        for subfam in subfamnames:
            sfindi=subfamnames.index(subfam)
            eval('node.add_features('+subfam+'_birth_rate='+str(allarray[sfindi][3][nodeindi])+')')
            eval('node.add_features('+subfam+'_death_rate='+str(allarray[sfindi][4][nodeindi])+')')
    annames=[]
    for subfam in subfamnames:
        annames.append(subfam+'_birth_rate')
        annames.append(subfam+'_death_rate')
    out=open(sys.argv[1],'w')
    out.write(write2nexus(antree,annames).replace('=0,','=0.0,').replace('=0]','=0.0]'))
    out.close()






        



