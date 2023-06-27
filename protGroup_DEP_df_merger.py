import pandas as pd

def file_merger(prot_group, dep_result):
    """Takes a proteinGroups df file and a DEP results df file and merges them on Gene name column"""

    # change the 'Gene name' of proteinGroups to 'name' for merge
    prot_group.rename(columns={"Gene names": "name"}, inplace=True)

    # update the 'name' column in proteinGroups to add .1, .2 etc for duplicate names (DEP does this)
    # will use a function that takes names, if it appears twice gets a .1, .2 etc

    def update_duplicate_names(names):
        name_count = {}
        updated_names = []
    
        for name in names:
            if name in name_count:
                name_count[name] += 1
                updated_name = f"{name}.{name_count[name]}"
            else:
                name_count[name] = 0
                updated_name = name
        
            updated_names.append(updated_name)
    
        return updated_names
    
    # Extract the 'names' column of proteinGroups as a list
    names = prot_group['name'].tolist()

    # Update duplicate names
    updated_name = update_duplicate_names(names)

    # Replace the 'names' column with the updated names
    prot_group['name'] = updated_name

    # merge the two results file together using name
    df_merged = pd.merge(prot_group, dep_result, how = 'inner', on = 'name')

    # write this to a csv file in the current directory
    df_merged.to_csv('merged_results_peptides.csv', index = False)


    
if __name__ == '__main__':

    # load in the user's proteinGroups.txt file and results.csv files as pandas df
    df_prot_group = pd.read_table(input("Please enter the proteinGroups.txt file: "))
    df_DEP = pd.read_csv(input("Please enter the csv results file: "), sep = ',', index_col=[0])

    # call the merger function on the two files
    file_merger(df_prot_group, df_DEP)

    






    
    
    



