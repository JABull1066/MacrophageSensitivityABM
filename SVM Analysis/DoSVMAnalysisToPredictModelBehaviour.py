import numpy as np
import seaborn as sns
import os
import pandas as pd
import matplotlib.pyplot as plt
sns.set_style('white')
sns.set(font_scale=5)

# Set path to ABM simulation outputs
pathToSims = './SimulationOutputs/'
pathToPCFData = './Data/' # Folder containing wPCFs and PCFs in .csv form

        
#%%
files = os.listdir(pathToPCFData) # Get all PCFs/wPCFs
allData = {}
allParams = {}
for file in files:
    print(file)
    data = np.genfromtxt(pathToPCFData + file, delimiter=',')
    allData[file] = data.flatten()
    temp = file.split('.')[0]
    temp = temp.split('_')
    ID = temp[0]
    time = temp[2].split('-')[1]
    nParams = temp[1][-1]
    params = pathToSims + f'ID-{ID}_time-{time}_From{nParams}ParamSweep_Params.csv'
    df_params = pd.read_csv(params)
    key = f'ID-{ID}_time-{time}_From{nParams}ParamSweep'
    allParams[key] = df_params
    
# Specify shape of the wPCF data
shape = (191,101)

# Train: iterations 1-5 of all of 81 param vals
trainSets = {v:allParams[v] for v in allParams if '2ParamSweep' in v if float(allParams[v]['tgfThresholdForPhenotypeSwitch'])==0.5 if int(allParams[v]['iterationNumber'])<5}
# Test: iterations 6-10 of all of 81 param vals
testSets = {v:allParams[v] for v in allParams if '2ParamSweep' in v if float(allParams[v]['tgfThresholdForPhenotypeSwitch'])==0.5 if int(allParams[v]['iterationNumber'])>4}
# Predict: all 6 parameter sweep simulations
predictSets = {v:allParams[v] for v in allParams if '6ParamSweep' in v}

df_train = pd.DataFrame()
for v in trainSets:
    tempDf = trainSets[v]
    tempDf['Filename'] = v
    df_train = pd.concat([df_train, tempDf], ignore_index=True, sort=False)
df_test = pd.DataFrame()
for v in testSets:
    tempDf = testSets[v]
    tempDf['Filename'] = v
    df_test = pd.concat([df_test, tempDf], ignore_index=True, sort=False)

df_predict = pd.DataFrame()
for v in predictSets:
    tempDf = predictSets[v]
    tempDf['Filename'] = v
    df_predict = pd.concat([df_predict, tempDf], ignore_index=True, sort=False)


pcfTypes = ['Vessel-wPCF','Tumour-wPCF','VesselToTumour-PCF']
pcfs_train = {v:{} for v in pcfTypes}
pcfs_test = {v:{} for v in pcfTypes}
pcfs_predict = {v:{} for v in pcfTypes}

for v in range(len(df_train)):
    row = df_train.iloc[v]
    sweepNumber = row.Filename.split('_')[2][4]
    for celltype in pcfTypes:
        key = f'{row.ID}_ParamSweep{sweepNumber}_Time-{row.Timepoint}_{celltype}.csv'
        if key in allData:
            pcfs_train[celltype][row.Filename] = allData[key]
for v in range(len(df_test)):
    row = df_test.iloc[v]
    sweepNumber = row.Filename.split('_')[2][4]
    for celltype in pcfTypes:
        key = f'{row.ID}_ParamSweep{sweepNumber}_Time-{row.Timepoint}_{celltype}.csv'
        if key in allData:
            pcfs_test[celltype][row.Filename] = allData[key]
for v in range(len(df_predict)):
    row = df_predict.iloc[v]
    sweepNumber = row.Filename.split('_')[2][4]
    for celltype in pcfTypes:
        key = f'{row.ID}_ParamSweep{sweepNumber}_Time-{row.Timepoint}_{celltype}.csv'
        if key in allData:
            pcfs_predict[celltype][row.Filename] = allData[key]

# Assign labels to parameter sets
labels = []
for v in range(len(df_train)):
    row = df_train.iloc[v]
    if row.chi_macrophageToCSF == 0.5:
        label = 'Equilibrium'
    elif row.chi_macrophageToCSF == 1:
        if row.halfMaximalExtravasationCsf1Conc > 0.3:
            label = 'Equilibrium'
        else:
            label = 'Elimination'
    else:
        if row.halfMaximalExtravasationCsf1Conc < 0.4:
            label = 'Elimination'
        elif row.halfMaximalExtravasationCsf1Conc == 0.4:
            if row.chi_macrophageToCSF < 3:
                label = 'Elimination'
            else:
                label = 'Escape'
        else:
            label = 'Escape'
    labels.append(label)
df_train['Label'] = labels

colDict = {'Equilibrium':plt.cm.tab10(0),'Elimination':plt.cm.tab10(2),'Escape':plt.cm.tab10(1)}
cols = [colDict[v] for v in df_train.Label]

# Repeat for testing data
labels = []
for v in range(len(df_test)):
    row = df_test.iloc[v]
    if row.chi_macrophageToCSF == 0.5:
        label = 'Equilibrium'
    elif row.chi_macrophageToCSF == 1:
        if row.halfMaximalExtravasationCsf1Conc > 0.3:
            label = 'Equilibrium'
        else:
            label = 'Elimination'
    else:
        if row.halfMaximalExtravasationCsf1Conc < 0.4:
            label = 'Elimination'
        elif row.halfMaximalExtravasationCsf1Conc == 0.4:
            if row.chi_macrophageToCSF < 3:
                label = 'Elimination'
            else:
                label = 'Escape'
        else:
            label = 'Escape'
    labels.append(label)
df_test['True Label'] = labels


#%% Combine PCFs and do dimension reduction / classification
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors

def ravelAndPlotWPCF(pcf,shape,cmap=plt.cm.plasma,cbarLabel=None):
    # Aid for visualisation of an unravelled wPCF
    vmin=0
    vmax = 8
    one = 1/vmax
    nCols = 1000
    threshold = one*nCols

    colors1 = plt.cm.Greens(np.linspace(0, 1, round(threshold)))
    colors2 = cmap(np.linspace(0, 1, nCols - round(threshold)))

    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    cmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    
    #plt.imshow(pcf.reshape(shape).T,origin='lower',cmap=cmap)
    plt.imshow(pcf.reshape(shape).T,origin='lower',extent=[0,1,0,1],cmap=cmap,vmin=vmin,vmax=vmax)
    ax = plt.gca()
    ax.set_xlabel('Radius $(r)$')
    tickProps = [0,0.25,0.5,0.75,1]
    ax.set_xticks(tickProps)
    ax.set_xticklabels([v*20 for v in tickProps])
    ax.set_ylabel('Phenotype $(p)$')
    plt.colorbar(label=cbarLabel)
      

    

dictionary = {}
for v in pcfs_train['Vessel-wPCF']:
    print(v)
    label = df_train[df_train.Filename == v].iloc[0].Label
    pcf = pcfs_train['Vessel-wPCF'][v]
    dictionary[v] = {'label':label,'pcf_v':pcf}
    if v in pcfs_train['Tumour-wPCF']:
        dictionary[v]['pcf_t'] = pcfs_train['Tumour-wPCF'][v]
    if v in pcfs_train['VesselToTumour-PCF']:
        dictionary[v]['pcf_vTot'] = pcfs_train['VesselToTumour-PCF'][v]

# Now expand each PCF signature out into a single array
Y = np.ones(shape=(len(dictionary),np.prod(shape)*2 + shape[0]))
labs = []
for i, v in enumerate(dictionary):
    labs.append(dictionary[v]['label'])
    Y[i,:np.prod(shape)] = dictionary[v]['pcf_v']
    if 'pcf_t' in dictionary[v]:
        Y[i,np.prod(shape):-shape[0]] = dictionary[v]['pcf_t']
    if 'pcf_vTot' in dictionary[v]:
        Y[i,-shape[0]:] = dictionary[v]['pcf_vTot']

# Y is now a n by 38773 matrix, where each of the n rows is the vectorised PCF signature for a single point cloud

labs = np.asarray(labs)

# impute missing values for Y - replace NaNs with a 1 (as this is the `uninformative' value for a PCF)
from sklearn.impute import SimpleImputer
imp_mean = SimpleImputer(missing_values=np.nan, strategy='constant',fill_value=1)
imp_mean.fit(Y)
Y = imp_mean.transform(Y)


# Choose number of principal components to reduce to - will require at least nComp different images to define
nComp = 100
from sklearn.decomposition import PCA
pca = PCA(n_components=nComp,svd_solver="randomized")
pca.fit(Y)

    
mu = np.mean(Y, axis=0)
Xhat = np.dot(pca.transform(Y)[:,:nComp], pca.components_[:nComp,:])
Xhat += mu



fig = plt.figure(figsize=(18, 18))
ax = fig.add_subplot(projection='3d')

for l in colDict:
    meanPCF = np.mean(Y[labs == l,:],axis=0)
    p_subset = pca.transform(Y)[labs==l,:nComp]
    ax.scatter(p_subset[:,0], p_subset[:,1], p_subset[:,2],c=colDict[l],s=50)
    mean_PCAspace = np.mean(p_subset,axis=0)
    ax.scatter(mean_PCAspace[0], mean_PCAspace[1], mean_PCAspace[2],c=[0.3,0.3,0.3],marker='D',s=200)
    ax.view_init(elev=20., azim=-121)
    ax.set_xlim([-150,600])
    ax.set_ylim([-0,600])
    ax.set_zlim([-250,250])
    ax.set_xlabel('\n\nPC1')
    ax.set_ylabel('\n\nPC2')
    ax.set_zlabel('PC3  ')

    Xhat_ofMeanInPCA = np.dot(mean_PCAspace[:nComp], pca.components_[:nComp,:])
    Xhat_ofMeanInPCA += mu
    # plt.figure(figsize=(18,18))
    # ravelAndPlotWPCF(meanPCF[0:np.prod(shape)],shape,plt.cm.plasma)
    # plt.savefig(f'./Figures/Fig7_PCA/{l}_MeanVesselPCF.svg')
    # plt.close()
    
    plt.figure(figsize=(18,18))
    ravelAndPlotWPCF(Xhat_ofMeanInPCA[0:np.prod(shape)],shape,plt.cm.plasma,cbarLabel='$wPCF(r,p,B)$')
    plt.savefig(f'./{l}_MeanVesselPCF.png')
    plt.close()
    
    # plt.figure(figsize=(18,18))
    # ravelAndPlotWPCF(meanPCF[np.prod(shape):-shape[0]],shape,plt.cm.inferno)
    plt.figure(figsize=(18,18))
    ravelAndPlotWPCF(Xhat_ofMeanInPCA[np.prod(shape):-shape[0]],shape,plt.cm.inferno,cbarLabel='$wPCF(r,p,T)$')
    plt.savefig(f'./{l}_MeanTumourPCF.png')
    plt.close()
    
    plt.figure(figsize=(18,18))
    # plt.plot(meanPCF[-shape[0]:],c=colDict[l],linestyle=':')
    plt.plot(0.1*np.arange(len(Xhat_ofMeanInPCA[-shape[0]:])),Xhat_ofMeanInPCA[-shape[0]:],c=colDict[l])
    plt.ylim([0,2.5])
    plt.ylabel('$g(r)$')
    plt.xlabel('Radius $(r)$')
    plt.gca().axhline(1,linestyle=':',c='k')
    plt.savefig(f'./{l}_MeanTtoVPCF.png')
    plt.close()
plt.show()
plt.savefig(f'./Figures/Fig7_PCA/PCA-space.png')

#%% SVM time
from sklearn import svm
clf = svm.SVC(kernel="rbf")

nComp = 100 # Number of principal components to use to define dimensionality of space for SVM classifier

P = pca.transform(Y)[:,:nComp]
clf.fit(P, labs)

dictionary_test = {}
for v in pcfs_test['Vessel-wPCF']:
    label = df_test[df_test.Filename == v].iloc[0]['True Label']
    pcf = pcfs_test['Vessel-wPCF'][v]
    dictionary_test[v] = {'label':label,'pcf_v':pcf}
    if v in pcfs_test['Tumour-wPCF']:
        dictionary_test[v]['pcf_t'] = pcfs_test['Tumour-wPCF'][v]
    if v in pcfs_test['VesselToTumour-PCF']:
        dictionary_test[v]['pcf_vTot'] = pcfs_test['VesselToTumour-PCF'][v]
# Now expand those out into an array
Y_test = np.ones(shape=(len(dictionary_test),np.prod(shape)*2 + shape[0]))
labs_test = []
for i, v in enumerate(dictionary_test):
    labs_test.append(dictionary_test[v]['label'])
    Y_test[i,:np.prod(shape)] = dictionary_test[v]['pcf_v']
    if 'pcf_t' in dictionary_test[v]:
        Y_test[i,np.prod(shape):-shape[0]] = dictionary_test[v]['pcf_t']
    if 'pcf_vTot' in dictionary_test[v]:
        Y_test[i,-shape[0]:] = dictionary_test[v]['pcf_vTot']

labs_test = np.asarray(labs_test)

# impute missing values for Y
from sklearn.impute import SimpleImputer
imp_mean = SimpleImputer(missing_values=np.nan, strategy='constant',fill_value=1)
imp_mean.fit(Y_test)
Y_test = imp_mean.transform(Y_test)
P_test = pca.transform(Y_test)[:,:nComp]

predicted_labels = clf.predict(P_test)
print(sum(labs_test == predicted_labels)/len(labs_test))


#%% Repeat process for `predict' dataset
dictionary_predict = {}
for v in pcfs_predict['Vessel-wPCF']:
    print(v)
    pcf = pcfs_predict['Vessel-wPCF'][v]
    dictionary_predict[v] = {'pcf_v':pcf}
    if v in pcfs_predict['Tumour-wPCF']:
        dictionary_predict[v]['pcf_t'] = pcfs_predict['Tumour-wPCF'][v]
    if v in pcfs_predict['VesselToTumour-PCF']:
        dictionary_predict[v]['pcf_vTot'] = pcfs_predict['VesselToTumour-PCF'][v]
# Now expand those out into an array
Y_predict = np.ones(shape=(len(dictionary_predict),np.prod(shape)*2 + shape[0]))

for i, v in enumerate(dictionary_predict):
    Y_predict[i,:np.prod(shape)] = dictionary_predict[v]['pcf_v']
    if 'pcf_t' in dictionary_predict[v]:
        Y_predict[i,np.prod(shape):-shape[0]] = dictionary_predict[v]['pcf_t']
    if 'pcf_vTot' in dictionary_predict[v]:
        Y_predict[i,-shape[0]:] = dictionary_predict[v]['pcf_vTot']


# impute missing values for Y
from sklearn.impute import SimpleImputer
imp_mean = SimpleImputer(missing_values=np.nan, strategy='constant',fill_value=1)
imp_mean.fit(Y_predict)
Y_predict = imp_mean.transform(Y_predict)
P_predict = pca.transform(Y_predict)[:,:nComp]
