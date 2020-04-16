#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from pandas import ExcelWriter
from io import StringIO
import io
import os

#Prepares the table for merger, concatenating the name with condition
def prep_table(file, sheet_name=0):
    data = pd.read_excel(file, sheet_name=sheet_name, header=None).dropna(axis=1, how='all')
    output = io.StringIO(data.to_csv(header=None, index=None, sep='\t'))

    output.seek(0)
    line1 = output.readline()
    line2 = output.readline()

    results = zip(line1.split('\t')[0::2],line1.split('\t')[1::2])
    line12 = list(results)
    results2 = zip(line2.split('\t')[0::2],line2.split('\t')[1::2])
    line22 = list(results2)

    final_line = []

    if len(line12) == len(line22):
        for x in range(len(line12)):
            first = line22[x][0]
            second = '##'.join([line12[x][0],line12[x][1]])
            final_line.append('\t'.join([first,second]))

    final_file = io.StringIO()
    output.seek(0)
    output.readline()
    output.readline()
    final_file.write('\t'.join(final_line))
    for line in output.readlines():
        final_file.write(line)
    final_file.seek(0)

    return pd.read_table(final_file)

#Merge dataframe with specified error in ppm
def merge(df, error=1):

    '''
    Merge and align m/z data from an prep_table(file) dataframe.

    Parameters:
        df : dataframe
        error : Error of aligment in ppm
    '''

    t = pd.DataFrame()
    for x in range(0,len(list(df.columns.values)),2):
        tempdf = df.iloc[:, [x,x+1]].dropna()
        tempdf.insert(2, 'Sample', list(df)[x+1])
        tempdf.columns = ['Mass', 'Intensity', 'Sample']
        t = pd.concat([t,tempdf])
    t = t.sort_values('Mass').reset_index(drop=True)
    t.Sample = t.Sample.astype('category')

    dfp = pd.DataFrame(columns=list(t.Sample.cat.categories))
    masses = []
    d = {}

    #Set first sample
    last = t.Mass[0]
    d[t.Sample[0]] = t.Intensity[0]
    masses.append(t.Mass[0])

    for row in t[1:].itertuples():
        if row[3] in d:
            c = round(sum(masses)/len(masses),6)
            dfp.loc[c] = pd.Series(d)
            masses = []
            d = {}
            last=row.Mass
            d[row[3]] = row[2]
            masses.append(row[1])

        elif abs((last-row[1])*1000000)/last <= error:
            d[row[3]] = row[2]
            masses.append(row[1])
        else:
            c = round(sum(masses)/len(masses),6)
            dfp.loc[c] = pd.Series(d)
            masses = []
            d = {}
            last=row[1]
            d[row[3]] = row[2]
            masses.append(row[1])

    if len(masses)==1:
        dfp.loc[last] = pd.Series(d)
    else:
        c = round(sum(masses)/len(masses),6)
        dfp.loc[c] = pd.Series(d)

    return dfp

#Split the dataframe based on the condition
def split_sample_label(df):
    ficheiro = io.StringIO(df.to_csv(sep='\t'))
    ficheiro.seek(0)
    first_liner = ficheiro.readline().strip('\n')
    work = first_liner.split('\t')[1:]

    line1 = ['Sample']
    line2 = ['Label']
    for x in work:
        line1.append(x.split('##')[0])
        line2.append(x.split('##')[1])

    new_file = io.StringIO()
    new_file.write('\t'.join(line2)+'\n')
    new_file.write('\t'.join(line1)+'\n')
    for line in ficheiro.readlines():
        new_file.write(line)
    new_file.seek(0)
    return(pd.read_table(new_file, header=[0,1], index_col=0))

# Clean the dataframe of low abundant metabolites. Masses must be present in at least "minsamples" (per condition) to be analysed
def cleandf(dataframe, minsamples=1):
    '''Clean a merged Dataframe based on the minimum samples'''

    labels = list(dataframe.columns.levels[0])
    labels

    for label in labels:
        dataframe['Counts',label] = dataframe[label].notnull().sum(axis=1)>=minsamples

    dataframe = dataframe.drop(dataframe[dataframe['Counts'].sum(axis=1)==0].index)
    dataframe.drop('Counts', axis=1, inplace=True)

    return dataframe.swaplevel(axis=1)
