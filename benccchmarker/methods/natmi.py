from functools import partial
import glob
import os
from pathlib import Path
import multiprocessing
import sys
import time
import tracemalloc
import yaml

import anndata
import pandas as pd
import scanpy as sc
import pandas as pd

def TransferToGeneSymbol(homoMapDir, speciestype, interSpeciestype, taxidCol, geneSymbolCol, hidCol, lrM):
    # load data
    homoDF = pd.read_csv(homoMapDir, sep='\t', index_col=None, header=None)
    # reduce the list to genes of the given species
    humanDF = homoDF.loc[homoDF[taxidCol] == int(interSpeciestype),]
    homoDF = homoDF.loc[homoDF[taxidCol] == int(speciestype),]

    ligandGIDList = list(lrM.index.values)
    receptorGIDList = list(lrM.columns.values)
    
    lhumanDF = humanDF.loc[humanDF[geneSymbolCol].isin(ligandGIDList), [hidCol, geneSymbolCol]]
    lhumanDF.columns = ['hid', 'gene']
    lhumanDF = lhumanDF.drop_duplicates(subset=['hid']).drop_duplicates(subset=['gene'])
    rhumanDF = humanDF.loc[humanDF[geneSymbolCol].isin(receptorGIDList), [hidCol, geneSymbolCol]]
    rhumanDF.columns = ['hid', 'gene']
    rhumanDF = rhumanDF.drop_duplicates(subset=['hid']).drop_duplicates(subset=['gene'])
    lrM = lrM.loc[lhumanDF.loc[:, 'gene'], rhumanDF.loc[:, 'gene']]
    lrM.index = lhumanDF.loc[:, 'hid']
    lrM.columns = rhumanDF.loc[:, 'hid']
    
    lhomoDF = homoDF.loc[homoDF[hidCol].isin(lhumanDF.loc[:, 'hid']), [hidCol, geneSymbolCol]]
    lhomoDF.columns = ['hid', 'gene']
    lhomoDF = lhomoDF.drop_duplicates(subset=['hid']).drop_duplicates(subset=['gene'])
    rhomoDF = homoDF.loc[homoDF[hidCol].isin(rhumanDF.loc[:, 'hid']), [hidCol, geneSymbolCol]]
    rhomoDF.columns = ['hid', 'gene']
    rhomoDF = rhomoDF.drop_duplicates(subset=['hid']).drop_duplicates(subset=['gene'])
    
    lrM = lrM.loc[lhomoDF.loc[:, 'hid'], rhomoDF.loc[:, 'hid']]
    lrM.index = lhomoDF.loc[:, 'gene']
    lrM.columns = rhomoDF.loc[:, 'gene']

    return lrM

def ClusterAnnotateEM(resultDir, emDF, ann):
    # get cluster list
    clusterIdList = sorted(list(set(ann.loc[:, 'cluster'].tolist())))
    
    # calculate the expressions for each cluster
    sumCounttableDFList = []
    meanCounttableDFList = []
    celltableDFList = []
    counttableDFList = []
    for clusterId in clusterIdList:
        # get the sub dataframe of the cluster
        cellsInClusterList = list(ann.index[ann['cluster'] == clusterId])
        clusterDF = emDF.loc[:, cellsInClusterList]
        
        # replace headers for the cluster
        sumDF = clusterDF.sum(axis=1).to_frame(name=clusterId)
        meanDF = clusterDF.mean(axis=1).to_frame(name=clusterId)
        # calculate the number of expressed cells
        cellDF = clusterDF[clusterDF>0].count(axis=1).astype(float).to_frame(name=clusterId)
        countDF = cellDF/len(cellsInClusterList)
        
        # add to the final dataframes
        sumCounttableDFList.append(sumDF)
        meanCounttableDFList.append(meanDF)
        celltableDFList.append(cellDF)
        counttableDFList.append(countDF)

    # merge results and save
    sumCounttableDF = pd.concat(sumCounttableDFList, axis=1)
    meanCounttableDF = pd.concat(meanCounttableDFList, axis=1)
    counttableDF = pd.concat(counttableDFList, axis=1)
    celltableDF = pd.concat(celltableDFList, axis=1)
    
    return sumCounttableDF, meanCounttableDF, counttableDF, celltableDF

def GenLigandReceptorList(pairsDF):
    ligandList = []
    receptorList = []
    ligandNameList = pairsDF.index.values
    for ligandName in ligandNameList:
        pairedRecptors = list(pairsDF.loc[ligandName, pairsDF.loc[ligandName,] > 0].index.values)
        subLigand = [ligandName] * len(pairedRecptors)
        ligandList = ligandList +subLigand
        receptorList = receptorList + pairedRecptors
    pairListDF = pd.DataFrame({'ligand':ligandList, 'receptor':receptorList})
    return pairListDF

def SplitIntoSinalProteins(sumEMDF, meanEMDF, countEMDF, cellEMDF, ligandApprovedSymbolList, receptorApprovedSymbolList):
    #split ligand and receptor
    ligandApprovedSymbolList = list(set(ligandApprovedSymbolList).intersection(set(sumEMDF.index)))
    receptorApprovedSymbolList = list(set(receptorApprovedSymbolList).intersection(set(sumEMDF.index)))
    sumligandDF = sumEMDF.loc[ligandApprovedSymbolList,:].fillna(0.0)
    sumreceptorDF = sumEMDF.loc[receptorApprovedSymbolList,:].fillna(0.0)
    meanligandDF = meanEMDF.loc[ligandApprovedSymbolList,:].fillna(0.0)
    meanreceptorDF = meanEMDF.loc[receptorApprovedSymbolList,:].fillna(0.0)
    countligandDF = countEMDF.loc[ligandApprovedSymbolList,:].fillna(0.0)
    countreceptorDF = countEMDF.loc[receptorApprovedSymbolList,:].fillna(0.0)
    cellligandDF = cellEMDF.loc[ligandApprovedSymbolList,:].fillna(0.0)
    cellreceptorDF = cellEMDF.loc[receptorApprovedSymbolList,:].fillna(0.0)
    
    # calculate specificity
    sumSpecifiedLigandDF = sumligandDF.div(sumligandDF.sum(axis=1), axis=0).fillna(0.0)
    sumSpecifiedReceptorDF = sumreceptorDF.div(sumreceptorDF.sum(axis=1), axis=0).fillna(0.0)
    meanSpecifiedLigandDF = meanligandDF.div(meanligandDF.sum(axis=1), axis=0).fillna(0.0)
    meanSpecifiedReceptorDF = meanreceptorDF.div(meanreceptorDF.sum(axis=1), axis=0).fillna(0.0)
    
    return cellligandDF, cellreceptorDF, countligandDF, countreceptorDF, sumligandDF, sumreceptorDF, meanligandDF, meanreceptorDF, sumSpecifiedLigandDF, sumSpecifiedReceptorDF, meanSpecifiedLigandDF, meanSpecifiedReceptorDF

def LRExpressions(typeString, ann, cellligandDF, cellreceptorDF, countligandDF, countreceptorDF, sumligandDF, sumreceptorDF, meanligandDF, meanreceptorDF, sumSpecifiedLigandDF, sumSpecifiedReceptorDF, meanSpecifiedLigandDF, meanSpecifiedReceptorDF, sourceFolder):
    origlabels = list(set(ann.iloc[:,0]))
    
    LRCountDict = {'Cluster':[],'Ligand count':[],'Receptor count':[]}
    
    writer = pd.ExcelWriter(os.path.join(sourceFolder, 'Ligands_Receptors_%s.xlsx' % (typeString.split('Edges_')[1])), engine='xlsxwriter')
    readmeDict = {'colA':['README'],'colB':['']}
    readmeDict['colA'].append('Ligand/receptor symbol')
    readmeDict['colB'].append('The official gene symbol of the detected ligand/receptor.')
    readmeDict['colA'].append('Total number of cells')
    readmeDict['colB'].append('Number of cells ligand/receptor is detected in.')
    readmeDict['colA'].append('Ligand/receptor detection rate')
    readmeDict['colB'].append('The ratio of cells that expressed the ligand/receptor to total cells in the cluster.')
    readmeDict['colA'].append('Ligand/receptor average expression value')
    readmeDict['colB'].append('The average expression level of the ligand/receptor in the cluster.')
    readmeDict['colA'].append('Ligand/receptor derived specificity of average expression value')
    readmeDict['colB'].append('The ratio of the average expression level of the ligand/receptor in the cluster to the sum of the average expression levels of the ligand/receptor in every cluster.')
    readmeDict['colA'].append('Ligand/receptor total expression value')
    readmeDict['colB'].append('The total expression level of the ligand/receptor in the cluster.')
    readmeDict['colA'].append('Ligand/receptor derived specificity of total expression value')
    readmeDict['colB'].append('The ratio of the total expression level of the ligand/receptor in the cluster to the sum of the total expression levels of the ligand/receptor in every cluster.')
    readmeDict['colA'].append('')
    readmeDict['colB'].append('')
    readmeDict['colA'].append('Summary of input file')
    readmeDict['colB'].append('')
    readmeDict['colA'].append('Cluster')
    readmeDict['colB'].append('Total number of cells')
    idmapDict = {}
    ididx = 0
    for origlabel in sorted(origlabels):
        ididx += 1
        idmapDict[origlabel] = 'Cluster %s' % ididx
        readmeDict['colA'].append(origlabel+'/'+idmapDict[origlabel])
        readmeDict['colB'].append(str(len(ann.index[ann['cluster'] == origlabel])))
    readmeDF = pd.DataFrame(readmeDict)
    readmeDF.to_excel(writer, sheet_name='README', index=False, header=False)
    
    for origlabel in sorted(origlabels):
        LRCountDict['Cluster'].append(origlabel)
        tempcellligandDF = cellligandDF.loc[:,[origlabel]]
        tempcountligandDF = countligandDF.loc[:,[origlabel]]
        tempsumligandDF = sumligandDF.loc[:,[origlabel]]
        tempmeanligandDF = meanligandDF.loc[:,[origlabel]]
        tempsumSpecifiedLigandDF = sumSpecifiedLigandDF.loc[:,[origlabel]]
        tempmeanSpecifiedLigandDF = meanSpecifiedLigandDF.loc[:,[origlabel]]
        tempLigandDF = pd.concat([tempcellligandDF, tempcountligandDF, tempmeanligandDF, tempmeanSpecifiedLigandDF, tempsumligandDF, tempsumSpecifiedLigandDF], axis=1)
        tempLigandDF = tempLigandDF.reset_index()
        tempLigandDF.columns = ['Ligand symbol', 'Total number of cells', 'Ligand detection rate', 'Ligand average expression value', 'Ligand derived specificity of average expression value', 'Ligand total expression value', 'Ligand derived specificity of total expression value']
        tempLigandDF['Total number of cells'] = tempLigandDF['Total number of cells'].astype(int)
        tempLigandDF = tempLigandDF.sort_values(by='Ligand detection rate', ascending=False)
        tempLigandDF = tempLigandDF.loc[tempLigandDF['Ligand total expression value']>0]
        LRCountDict['Ligand count'].append(len(tempLigandDF))
        tempcellreceptorDF = cellreceptorDF.loc[:,[origlabel]]
        tempcountreceptorDF = countreceptorDF.loc[:,[origlabel]]
        tempsumreceptorDF = sumreceptorDF.loc[:,[origlabel]]
        tempmeanreceptorDF = meanreceptorDF.loc[:,[origlabel]]
        tempmeanSpecifiedReceptorDF = meanSpecifiedReceptorDF.loc[:,[origlabel]]
        tempsumSpecifiedReceptorDF = sumSpecifiedReceptorDF.loc[:,[origlabel]]
        tempReceptorDF = pd.concat([tempcellreceptorDF, tempcountreceptorDF,tempmeanreceptorDF,tempmeanSpecifiedReceptorDF, tempsumreceptorDF,tempsumSpecifiedReceptorDF], axis=1)
        tempReceptorDF = tempReceptorDF.reset_index()
        tempReceptorDF.columns = ['Receptor symbol', 'Total number of cells', 'Receptor detection rate', 'Receptor average expression value', 'Receptor derived specificity of average expression value', 'Receptor total expression value', 'Receptor derived specificity of total expression value']
        tempReceptorDF['Total number of cells'] = tempReceptorDF['Total number of cells'].astype(int)
        tempReceptorDF = tempReceptorDF.sort_values(by='Receptor detection rate', ascending=False)
        tempReceptorDF = tempReceptorDF.loc[tempReceptorDF['Receptor total expression value']>0]
        LRCountDict['Receptor count'].append(len(tempReceptorDF))
        if len(tempLigandDF) > 0:
            tempLigandDF.to_excel(writer, sheet_name='Ligands in %s' % idmapDict[origlabel], index=False, header=True, columns = ['Ligand symbol', 'Total number of cells', 'Ligand detection rate', 'Ligand average expression value', 'Ligand derived specificity of average expression value', 'Ligand total expression value', 'Ligand derived specificity of total expression value'])
            worksheet = writer.sheets['Ligands in %s' % idmapDict[origlabel]]
            worksheet.conditional_format('C2:C%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('D2:D%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('E2:E%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('F2:F%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('G2:G%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            
        if len(tempReceptorDF) > 0:
            tempReceptorDF.to_excel(writer, sheet_name='Receptors in %s' % idmapDict[origlabel], index=False, header=True, columns = ['Receptor symbol', 'Total number of cells', 'Receptor detection rate', 'Receptor average expression value', 'Receptor derived specificity of average expression value', 'Receptor total expression value', 'Receptor derived specificity of total expression value'])
            worksheet = writer.sheets['Receptors in %s' % idmapDict[origlabel]]
            worksheet.conditional_format('C2:C%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('D2:D%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('E2:E%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('F2:F%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
            worksheet.conditional_format('G2:G%s' % (len(tempLigandDF)+1), {'type':'2_color_scale','min_color': '#FFFFFF','max_color': '#FF0000'})
    writer.save()
    

def FindCellsOfProtein(protein, proteinType, cellDF, countDF, sumDF, meanDF, sumSpecifiedDF, meanSpecifiedDF):
    subCellDF = cellDF.loc[cellDF.index == protein, :].T
    subCountDF = countDF.loc[countDF.index == protein, :].T
    subSumDF = sumDF.loc[sumDF.index == protein, :].T
    subMeanDF = meanDF.loc[meanDF.index == protein, :].T
    subSumSpecifiedDF = sumSpecifiedDF.loc[sumSpecifiedDF.index == protein, :].T
    subMeanSpecifiedDF = meanSpecifiedDF.loc[meanSpecifiedDF.index == protein, :].T
    mergedProteinDF = pd.concat([subCellDF, subCountDF, subSumDF, subMeanDF, subSumSpecifiedDF, subMeanSpecifiedDF], axis=1)
    newHeaders = ['cell '+proteinType, 'frequency '+proteinType, 'sum '+proteinType, 'mean '+proteinType, 'specified sum '+proteinType, 'specified mean '+proteinType]
    mergedProteinDF.columns = newHeaders
    
    return mergedProteinDF

def BuildHalfEdge(protein, cellDF, countDF, sumDF, meanDF, sumSpecifiedDF, meanSpecifiedDF, cellRole, proteinType):
    mergedProteinDF = FindCellsOfProtein(protein, proteinType, cellDF, countDF, sumDF, meanDF, sumSpecifiedDF, meanSpecifiedDF)
    
    #update headers
    mergedProteinDF[cellRole+' name'] = mergedProteinDF.index.values
    mergedProteinDF[proteinType] = protein
    mergedProteinDF.reset_index(drop=True, inplace=True)
    
    return mergedProteinDF

def GenSingleCell2CellEdge(ligandApprovedSymbolDict, receptorApprovedSymbolDict, resultDir, pair):
    ligand = pair[0]
    receptor = pair[1]
    if ligand not in ligandApprovedSymbolDict.keys() or receptor not in receptorApprovedSymbolDict.keys():
        return
    
    fn = os.path.join(resultDir,'%s-%s.xlsx' % (pair[0],pair[1]))
    
    edgeList = []
    ligandCellListDF = ligandApprovedSymbolDict[ligand]
    receptorCellListDF = receptorApprovedSymbolDict[receptor]
    for ligandCellListIndex in ligandCellListDF.index:
        for receptorCellListIndex in receptorCellListDF.index:
            ligandCellDF = ligandCellListDF.loc[ligandCellListIndex, :].to_frame().T.reset_index(drop=True)
            receptorCellDF = receptorCellListDF.loc[receptorCellListIndex, :].to_frame().T.reset_index(drop=True)
            mergedPairDF = pd.concat([ligandCellDF, receptorCellDF], axis=1)
            edgeList.append(mergedPairDF)
    
    edgeListDF = pd.concat(edgeList, axis=0)
    #add prodect columns
    edgeListDF['product of sum'] = edgeListDF['sum ligand'] * edgeListDF['sum receptor']
    edgeListDF['product of mean'] = edgeListDF['mean ligand'] * edgeListDF['mean receptor']
    edgeListDF['product of specified sum'] = edgeListDF['specified sum ligand'] * edgeListDF['specified sum receptor']
    edgeListDF['product of specified mean'] = edgeListDF['specified mean ligand'] * edgeListDF['specified mean receptor']
    
    #remove values equal zero
    edgeListDF = edgeListDF[edgeListDF['product of sum'] > 0]
    if len(edgeListDF) > 0:
        edgeListDF = edgeListDF.loc[:,['sending cluster name', 'ligand', 'receptor', 'target cluster name', 'cell ligand', 'frequency ligand', 'mean ligand', 'sum ligand', 'specified mean ligand', 'specified sum ligand', 
           'cell receptor', 'frequency receptor', 'mean receptor', 'sum receptor', 'specified mean receptor', 'specified sum receptor', 'product of mean', 'product of sum', 'product of specified mean', 'product of specified sum']]
        newCols = ['Sending cluster', 'Ligand symbol', 'Receptor symbol', 'Target cluster', 'Ligand-expressing cells', 'Ligand detection rate', 'Ligand average expression value', 'Ligand total expression value', 
            'Ligand derived specificity of average expression value', 'Ligand derived specificity of total expression value',  'Receptor-expressing cells', 'Receptor detection rate', 'Receptor average expression value', 'Receptor total expression value', 
            'Receptor derived specificity of average expression value', 'Receptor derived specificity of total expression value', 'Edge average expression weight', 'Edge total expression weight', 
            'Edge average expression derived specificity', 'Edge total expression derived specificity']
        edgeListDF.columns = newCols
        edgeListDF.to_excel(fn, header=True,index=False,columns=newCols)
    return

def GenCell2CellEdges(interDB, pairListDF, ligandApprovedSymbolDict, receptorApprovedSymbolDict, coreNum, resultDir):
    pairList = [tuple(x) for x in pairListDF.values]
    fn = os.path.join(resultDir,'LR-pairs_'+interDB)
    if not os.path.exists(fn):
        os.mkdir(fn)
    
    p = multiprocessing.Pool(coreNum)
    func = partial(GenSingleCell2CellEdge, ligandApprovedSymbolDict, receptorApprovedSymbolDict, fn)
    p.map(func, pairList)
    p.close()
    p.join()
    
    fl = []
    for p in pairList:
        f = os.path.join(fn,'%s-%s.xlsx' % (p[0],p[1]))
        if not os.path.exists(f):
            continue
        d = pd.read_excel(f, index_col=False, header=0)
        fl.append(d)
    edgeListDF = pd.concat(fl, axis=0)
    print('#### collected %s LR-mediated edges' % len(edgeListDF))
    
    return edgeListDF

def GenerateCell2CellTable(interDB, sumEMDF, meanEMDF, countEMDF, cellEMDF, typeString, pairsDF, ann, resultDir, coreNum):
    # generate ligand-receptor pair list
    pairListDF = GenLigandReceptorList(pairsDF)
    ligandApprovedSymbolList = list(pairsDF.index.values)
    receptorApprovedSymbolList = list(pairsDF.columns.values)
    
    print('#### extract signaling factors')
    # split original and speciefied countable to ligand and receptor file
    cellligandDF, cellreceptorDF, countligandDF, countreceptorDF, sumligandDF, sumreceptorDF, meanligandDF, meanreceptorDF, sumSpecifiedLigandDF, sumSpecifiedReceptorDF, meanSpecifiedLigandDF, meanSpecifiedReceptorDF = SplitIntoSinalProteins(sumEMDF, meanEMDF, countEMDF, cellEMDF, ligandApprovedSymbolList, receptorApprovedSymbolList)
    LRExpressions(typeString, ann, cellligandDF, cellreceptorDF, countligandDF, countreceptorDF, sumligandDF, sumreceptorDF, meanligandDF, meanreceptorDF, sumSpecifiedLigandDF, sumSpecifiedReceptorDF, meanSpecifiedLigandDF, meanSpecifiedReceptorDF, resultDir)
    
    # find cells for ligand and receptor
    ligandApprovedSymbolDict = {}
    for ligandApprovedSymbol in sumligandDF.index:#ligandApprovedSymbolList:
        mergedLigandDF = BuildHalfEdge(ligandApprovedSymbol, cellligandDF, countligandDF, sumligandDF, meanligandDF, sumSpecifiedLigandDF, meanSpecifiedLigandDF, 'sending cluster', 'ligand')
        ligandApprovedSymbolDict[ligandApprovedSymbol] = mergedLigandDF
    receptorApprovedSymbolDict = {}
    for receptorApprovedSymbol in sumreceptorDF.index:#receptorApprovedSymbolList:
        mergedReceptorDF = BuildHalfEdge(receptorApprovedSymbol, cellreceptorDF, countreceptorDF, sumreceptorDF, meanreceptorDF, sumSpecifiedReceptorDF, meanSpecifiedReceptorDF, 'target cluster', 'receptor')
        receptorApprovedSymbolDict[receptorApprovedSymbol] = mergedReceptorDF
    print('#### construct signaling interactions')
    
    # construct and save cell-to-cell edge
    edgeListDF = GenCell2CellEdges(interDB, pairListDF, ligandApprovedSymbolDict, receptorApprovedSymbolDict, coreNum, resultDir)
    
    if len(edgeListDF) > 0:
        return edgeListDF
    else:
        return []
        print('#### No signaling interactions found')

def GenerateDataFiles(interDB, sumEMDF, meanEMDF, countEMDF, cellEMDF, typeString, lrM, ann, resultDir, coreNum):
    # generate cell-ligand-receptor-cell table
    edgeListDF = GenerateCell2CellTable(interDB, sumEMDF, meanEMDF, countEMDF, cellEMDF, typeString, lrM, ann, resultDir, coreNum)
    
    #save result
    if len(edgeListDF) > 0:
        saveann = ann.reset_index()
        saveann.columns = ['cell', 'cluster']
        saveann.to_csv(os.path.join(resultDir,'ClusterMapping.csv'), index=False, header=True)
    
        edgeListFileDir = os.path.join(resultDir, typeString + '.csv')
        edgeListDF.to_csv(edgeListFileDir, index=False, header=True,columns=['Sending cluster', 'Ligand symbol', 'Receptor symbol', 'Target cluster', 'Ligand-expressing cells', 'Ligand detection rate', 'Ligand average expression value', 'Ligand total expression value', 
                'Ligand derived specificity of average expression value', 'Ligand derived specificity of total expression value',  'Receptor-expressing cells', 'Receptor detection rate', 'Receptor average expression value', 'Receptor total expression value', 
                'Receptor derived specificity of average expression value', 'Receptor derived specificity of total expression value', 'Edge average expression weight', 
                'Edge average expression derived specificity', 'Edge total expression weight', 'Edge total expression derived specificity'])
        with open(os.path.join(resultDir,'README.txt'), 'w') as file_object:
            file_object.write('README\n')
            file_object.write('\n')
            file_object.write('ClusterMapping.csv: cell-to-cluster mapping.\n')
            file_object.write('Ligands_Receptors_xx.xlsx: information about ligands and receptors in each cell-type/single-cell cluster.\n')
            file_object.write('Edges_xx.csv: all ligand-receptor-mediated communications.\n')
            file_object.write('LR-pairs_xx folder: all ligand-receptor-mediated communications via a ligand-receptor pair.\n')
            file_object.write('\n')
            file_object.write('Sending cluster: the cluster that expresses the ligand.\n')
            file_object.write('Target cluster: the cluster that expresses the receptor.\n')
            file_object.write('Ligand/receptor symbol: the official gene symbol of the detected ligand/receptor.\n')
            file_object.write('Ligand/receptor-expressing cells: number of cells ligand/receptor is detected in.\n')
            file_object.write('Ligand/receptor detection rate: the ratio of cells that expressed the ligand/receptor to total cells in the cluster.\n')
            file_object.write('Ligand/receptor average expression value: the average expression level of the ligand/receptor in the cluster.\n')
            file_object.write('Ligand/receptor derived specificity of average expression value: the ratio of the average expression level of the ligand/receptor in the cluster to the sum of the average expression levels of the ligand/receptor in every cluster.\n')
            file_object.write('Ligand/receptor total expression value: the total expression level of the ligand/receptor in the cluster.\n')
            file_object.write('Ligand/receptor derived specificity of total expression value: the ratio of the total expression level of the ligand/receptor in the cluster to the sum of the total expression levels of the ligand/receptor in every cluster.\n')
            file_object.write('Edge average expression weight: the product of average expression levels of the ligand and the receptor in the corresponding cluster(s).\n')
            file_object.write('Edge average expression derived specificity: the product of the ligand and receptor derived specificity of average expression values.\n')
            file_object.write('Edge total expression weight: the product of total expression levels of the ligand and the receptor in the corresponding cluster(s).\n')
            file_object.write('Edge total expression derived specificity: the product of the ligand and receptor derived specificity of total expression values.\n')
    else:
        with open(os.path.join(resultDir,'README.txt'), 'w') as file_object:
            file_object.write('README\n')
            file_object.write('\n')
            file_object.write('No signaling interactions found')

def main():
    # Get arguments from the command line
    args = sys.argv

    tracemalloc.start()
    start_time = time.time()

    current_file = Path(__file__)

    adata = anndata.read_h5ad(args[1])

    config = yaml.safe_load(open(args[3]))

    # Normalise data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    natmidb_file_path = current_file.parent.parent / "datasets" / "defaults" / "natmi" / "database" / "lrc2p.csv"

    main_output_path = "/".join(args[4].split("/")[:-1])
    output_path = f"{main_output_path}/tmp/natmi"

    os.makedirs(output_path, exist_ok=True)
    
    counts_file_path = f"{output_path}/em.csv"
    meta_file_path = f"{output_path}/metadata.csv"

    adata.to_df().T.to_csv(counts_file_path, index=True, header=True)
    adata.obs['cell_type'].to_csv(meta_file_path, index=True, header=True)

    em = pd.read_csv(counts_file_path, index_col=0, header=0)
    ann = pd.read_csv(meta_file_path, index_col=0, header=0)

    ann.columns = ['cluster']
    if len(set(em.columns).intersection(set(ann.index))) == 0:
        if len(set(em.columns).intersection(set(ann['cluster']))) == 0:
            raise ValueError("Cannot find matched annotations for the barcodes.")
        elif len(set(em.columns).intersection(set(ann['cluster']))) < len(set(em.columns)):
            raise ValueError("Cannot find matched annotations for all barcodes.")
        else:
            ann.columns = ['new.barcode']
            ann.reset_index().set_index('new.barcode')
            ann.columns = ['cluster']
    elif len(set(em.columns).intersection(set(ann.index))) < len(set(em.columns)):
        raise ValueError("Cannot find matched annotations for all barcodes.")
    
    lrL = pd.read_csv(natmidb_file_path, dtype=object, index_col=None, header=0)

    lset = sorted(list(set(lrL['Ligand gene symbol'])))
    rset = sorted(list(set(lrL['Receptor gene symbol'])))
    lrM = pd.DataFrame(0,index=lset,columns=rset)
    for idx in lrL.index:
        lrM.loc[lrL.loc[idx,'Ligand gene symbol'], lrL.loc[idx,'Receptor gene symbol']] = 1

    end_preprocessing_time = time.time()

    coreNum = 1
    sumEMDF, meanEMDF, countEMDF, cellEMDF = ClusterAnnotateEM(output_path, em, ann)

    GenerateDataFiles("lr2cp", sumEMDF, meanEMDF, countEMDF, cellEMDF, 'Edges_%s' % (os.path.basename("lrc2p")), lrM, ann, output_path, coreNum)

    end_inference_time = time.time()

    result_df = pd.read_csv(f"{output_path}/Edges_lrc2p.csv")

    result_df.rename(columns={"Sending cluster": "source", "Target cluster": "target", "Ligand symbol": "ligand", "Receptor symbol": "receptor"}, inplace=True)

    result_df["label"] = result_df["source"] + "---" + result_df["ligand"] + "---" + result_df["target"] + "---" + result_df["receptor"]

    result_df.to_csv(args[4], index=False)
            
    _, peak = tracemalloc.get_traced_memory()

    met_path = args[4].replace(".csv", "_met.csv")

    # Create a dataframe to store the memory usage in a column called peak_ram_used_mib
    preprocessing_time = end_preprocessing_time - start_time
    inference_time = end_inference_time - end_preprocessing_time

    met_df = pd.DataFrame(data={"peak_ram_used_mib": [peak], "preprocessing_time": [preprocessing_time], "inference_time": [inference_time]})
    met_df.to_csv(met_path, index=False)

if __name__ == "__main__":
    main()

