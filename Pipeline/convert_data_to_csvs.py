import os
header = "Index,Time Interval,Start,Died,Lost to Follow-up,Observed Interval,Observed Cum,Expected Interval,Expected Cum,Relative Interval,Relative Cum,SE Obs Interval,SE Obs Cum,SE Relative Interval,SE Relative Cum"
name_of_csv_relative_i = -2
i = 2
out_csvs = []
with open(snakemake.input[0]) as h:
    all_data = h.readlines()


while i < len(all_data):
    cur_title = all_data[i + name_of_csv_relative_i].replace("/","~").replace(",","").replace("\n", "")
    
    cur_data_list = []
    
    
    while i < len(all_data):
        i += 1
        if all_data[i].split(",")[0].isdigit():
            cur_data_list.append(all_data[i])
        else:
            break
    cur_data = "".join(cur_data_list)
    
    out_csvs.append((cur_title, header+"\n"+cur_data))
    
    while True:
        i+=1
        if i >= len(all_data) or (",,,,,,,,,,,,,," in all_data[i] and len(all_data[i]) == 15):
            break
    i+=3

for cur_title, cur_data in out_csvs:
    if os.path.isdir(snakemake.output[0]):
        os.mkdir(snakemake.output[0])
    with open(f"{snakemake.output[0]}/{cur_title}.csv", "w") as h:
        h.write(cur_data)
