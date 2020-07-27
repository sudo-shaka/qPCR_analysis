#!/usr/bin/python
try:
    import pandas as pd
    import statistics
    import numpy as np
    import matplotlib.pyplot as plt
except:
    print("Missing Required Dependancies")
    

csv_file = input("Enter the path to CSV file: ")

try:
    values = pd.read_csv(csv_file)

    data = pd.DataFrame(values,columns=['Sample Name','Target Name','CT'])

    samples = data.drop_duplicates(['Sample Name'])
    sample_name = samples['Sample Name'].values.tolist()
    gene_name = data.drop_duplicates(['Target Name'])
    target_names = gene_name['Target Name'].values.tolist()

    ave_data = data.drop_duplicates(['Sample Name','Target Name'],keep= 'last').reset_index(drop=True)
    list_data = data['CT'].values.tolist(); list_ave_data = ave_data['CT'].values.tolist()

except:
    print("Invalid file name..  Please try again!"); quit()

def AVE_CT():
#average CT values
    value = []; n_reps = len(list_data)/len(list_ave_data)
    try:
        for i in range(0,len(list_data)):
            for j in range(0,len(target_names)):
                for k in range(0,len(sample_name)):
                    if str(data.iloc[i,[0]].values) == str(samples.iloc[k,[0]].values) and str(data.iloc[i,[1]].values) == str(gene_name.iloc[j,[1]].values):
                        val = float(data.iloc[i,[2]])
                        value.append(val)
                        ave_data.loc[(j+(len(target_names))*(k-1)+len(target_names)),'CT'] = statistics.mean(value)
                    elif len(value) >= n_reps:
                        value=[]
    except:
        print("[!!] Error averging CT values"); quit()

#getting controls
def get_controls():
    global control_sample, control_target
    control_sample = sample_name[0]; control_target = target_names[0]
    print("\nUsing", control_sample,"and",control_target, "as normalization controls\n")

    correct = input("Are these parameters correct? (y/n) >> ").lower()
    if correct == "n":
        for sample in range(0,len(sample_name)):
            print('{:^3s}{:10s}'.format(str(sample),sample_name[sample]))
        try:
            n = int(input("What number corresponds to your control sample? >> "))
            control_sample = sample_name[int(n)]
        except:
            print("Please enter an integer...")
            get_controls()

        for control in range(0,len(target_names)):
            print('{:^3s}{:10s}'.format(str(control),target_names[control]))
        try:
            n = int(input("What number corresponds to your control gene? >> "))
            control_target = target_names[int(n)]
        except:
            print("Please enter an integer")
            get_controls()
        
        print("\nNow using", control_sample, "and", control_target, "as normalization controls... ")


#Getting DCT values
def get_DCT():
    try:
        for i in range(0,len(list_ave_data)):
                if ave_data.loc[i,'Target Name'] == control_target:
                    ave_cntl_ct = ave_data.loc[i,'CT']
                    ave_data.loc[i,'DCT'] = ""
                else:
                    ave_data.loc[i,'DCT'] = ave_data.loc[i,'CT'] - ave_cntl_ct
    except:
        print("[!!] Error caluclating DCT values")

#Getting DDCT values
def get_DDCT():
    DCT = [None]*len(target_names)
    for i in range(0,len(list_ave_data)):
        if ave_data.loc[i,'DCT'] != "":
            for j in range(0,len(target_names)):
                if target_names[j] == ave_data.loc[i,'Target Name']:
                    if ave_data.loc[i,'Sample Name'] == control_sample:
                        DCT[j] = ave_data.loc[i,'DCT']
                    else:
                        ave_data.loc[i,'DDCT'] = ave_data.loc[i,'DCT'] - DCT[j]
        else:
            ave_data.loc[i,'DDCT'] = ""

#Getting fold change
def get_Fold_Change():
    try:
        for i in range(0,len(list_ave_data)):
            if ave_data.loc[i,'DDCT'] == "":
                ave_data.loc[i,'Fold Change'] = ""
            elif str(ave_data.loc[i,'DDCT']) == 'nan':
                ave_data.loc[i,'Fold Change'] = 1
            else:
                ave_data.loc[i,'Fold Change'] = 2**(-1*ave_data.loc[i,'DDCT'])
    except:
        print("[!!] Error getting Fold Change")

    #Determining Change
    try:
        for i in range(0,len(list_ave_data)):
            try:
                if float(ave_data.loc[i,'Fold Change']) > 1.5:
                    ave_data.loc[i,'Up/Down'] = "Up"
                elif float(ave_data.loc[i,'Fold Change']) < 0.5:
                    ave_data.loc[i,'Up/Down'] = "Down"
                else:
                    ave_data.loc[i,'Up/Down'] = None
            except:
                ave_data.loc[i,'Up/Down'] = ""
    except:
        print("[!] Error : Cannot determine Fold Change Information")

def save_data():
    try:
        print(ave_data)
        csv_name = str(input('\nEnter name/path you would like to save data as: '))
        ave_data.to_csv(csv_name, encoding='utf-8', index=False)
        print("\nSaved Data!")
    except:
        print("\n[!!] Error: Unable to save CSV")

def createList(r1, r2):
    return np.arange(r1, r2+1, 1)

def plotting():
    hieght = []; xval = []; target = []; targ = ave_data.loc[0,'Target Name']; left = []
    try:
        for j in range(0,len(target_names)):
            for i in range(1,len(list_ave_data)):
                if ave_data.loc[i,'Fold Change'] != "" and ave_data.loc[i,'Target Name'] == target_names[j]:
                    ex = ave_data.loc[i,'Fold Change']
                    hieght.append(ex)
                    xv = ave_data.loc[i,'Sample Name']
                    xval.append(xv)
                    targ_b = targ
                    targ = ave_data.loc[i,'Target Name']
                    target.append(targ)
                    if targ != targ_b and len(hieght) > 1:
                        r1 = 1; r2 = len(hieght); left = createList(r1,r2); tick_label = xval
                        plt.bar(left,hieght,tick_label=tick_label, width = 0.8,color = 'black')
                        plt.ylabel(targ_b); plt.show()
                        hieght = [] ; xval = []
        r1 = 1; r2 = len(hieght); left = createList(r1,r2); tick_label = xval
        plt.bar(left,hieght,tick_label = tick_label, width =0.8, color = 'black')
        plt.ylabel(targ_b); plt.show()
    except:
        print("[!!] Error plotting")

AVE_CT()
get_controls()
get_DCT()
get_DDCT()
get_Fold_Change()
print(ave_data)


a_save = input("Would you like to save this data to a CSV file? >> ").lower()
if a_save == 'y' or a_save == 'yes':
    save_data()

a_plot = input("Would you like to plot this data? (Y/n) : ").lower()
if a_plot == 'y' or a_plot == 'yes':
    plotting()

    
