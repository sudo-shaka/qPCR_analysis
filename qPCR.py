#!/usr/bin/python3
try:
    import numpy as np
    import pandas as pd
    import argparse
except:
    print("Missing required dependances")
    quit()

def ERROR():
    print("[!]"); quit()

def IS_IN(len,VALUE,LIST):
    found = False
    for i in LIST :
        if i == VALUE:
            found = True
            break
    return found

def GET_DATA():
    ap = argparse.ArgumentParser()
    ap.add_argument("-d","--data",required = True, help = "Path to data")
    args = vars(ap.parse_args())
    filename = args['data']

    DATA = pd.read_csv(filename)
    DATA = pd.DataFrame(DATA, columns = ['Sample Name','Target Name','CT'])

    return DATA

def GET_MEANS(DATA):
    CONDITIONS = []; n_cond = 0
    for i in range(0,len(DATA)):
        CONDITIONS.append(str(DATA.iloc[i,0:2].values))
    
    for i in range(0,len(DATA)):
        if  i==0 or CONDITIONS[i] not in CONDITIONS[:i]:
            DATA.loc[i,'AVE_CT'] = DATA.loc[i,'CT']
            DATA.loc[i,'n_rep'] = 1; n_cond += 1
        else:
            for j in range(0,len(DATA)):
                if CONDITIONS[j] == CONDITIONS[i]:
                    DATA.loc[j,'n_rep'] += 1
                    if (DATA.loc[i,'CT'] - DATA.loc[j,'AVE_CT'] < -1 or\
                        DATA.loc[i,'CT'] - DATA.loc[j,'AVE_CT'] > 1):
                        print("Warning: {} condition's CT values have a differnce > 1. Could mean internal replicates are inconsistent.".format(CONDITIONS[i]))
                    DATA.loc[j,'AVE_CT'] = (DATA.loc[j,'AVE_CT']*DATA.loc[j,'n_rep'] + DATA.loc[i,'CT']) /  (DATA.loc[j,'n_rep'] + 1)
                    break
    DATA = DATA.dropna().reset_index(drop = True)
    
    return DATA

def GET_DCT(i,control_location,DATA):
    if DATA.loc[i,'Sample Name'] == DATA.loc[control_location,'Sample Name']:
        if DATA.loc[i,'Target Name'] == DATA.loc[control_location,'Target Name']:
            DCT = 0
        else:
            DCT = DATA.loc[i,'CT'] - DATA.loc[control_location,'CT']
    else:
        for j in range(0,len(DATA)):
            if DATA.loc[i,'Sample Name'] == DATA.loc[j,'Sample Name'] and \
                DATA.loc[j,'Target Name'] == DATA.loc[control_location, 'Target Name']:
                DCT = DATA.loc[i,'CT'] - DATA.loc[j,'CT']
                break
                
    return DCT
def GET_DDCT(i, control_location, DATA):
    DDCT = None
    if DATA.loc[i,'Sample Name'] == DATA.loc[control_location,'Sample Name'] or \
        DATA.loc[i,'Target Name'] == DATA.loc[control_location,'Target Name']:
        DDCT = 0
    else:
        for j in range(0,len(DATA)):
            if(DATA.loc[j,'Target Name']) == DATA.loc[i,'Target Name'] and \
                DATA.loc[j,'Sample Name'] == DATA.loc[control_location,'Sample Name']:
                DDCT = DATA.loc[i,'DCT'] - DATA.loc[j,'DCT']
                break
    return DDCT

def main():
    DATA = GET_DATA()
    DATA = GET_MEANS(DATA)
    control_location=0
    ask_correct = str(input("Is {} and {} your control condition? (y/n) >> ".format(DATA.loc[0,'Sample Name'], DATA.loc[0,'Target Name'])))
    if ask_correct[0].lower() != "y":
        print(DATA)
        control_location = int(input("Enter the number that corresponds to the control conditon >> "))

    print("Calculating... ")

    for i in range(0,len(DATA)):
        DATA.loc[i,'DCT'] = GET_DCT(i,control_location,DATA)
        DATA.loc[i,'DDCT'] = GET_DDCT(i,control_location,DATA)
        DATA.loc[i,'FOLD_CHANGE'] = 2**(-1*DATA.loc[i,'DDCT'])
    
    print(DATA)
    print("Done! Saving file to Results.csv")

    try:
        DATA.to_csv('Results.csv')
    except:
        ERROR()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        quit()