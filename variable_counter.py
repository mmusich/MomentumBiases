with open('passed_regions.txt', 'r') as f:
    lines = f.readlines() # read and store all lines into list, note list index starts from 0

lines_good = []
with open('passed_regions_edit.txt', 'w') as fout:
    for line in enumerate(lines):
        rev = line[1].split("_")[2] +"_"+ line[1].replace('\n', '').split("_")[3] +"_"+ line[1].split("_")[0] +"_"+ line[1].split("_")[1]
        flag = 0
        for line_tbc in enumerate(lines):
            if line_tbc[1].replace('\n', '') == rev:
                flag = 1 
                if line[0] < line_tbc[0]:
                    lines_good.append(line[1].replace('\n', ''))
                    print(line[1].replace('\n', ''), file=fout)
        if flag == 0:
            lines_good.append(line[1].replace('\n', ''))
            print(line[1].replace('\n', ''), file=fout)

#generate the IDs for all regions e.g 0_0, 0_1 etc. and check if there are any IDs that do not appear twice
nbinseta = 24
nbinspt = 5
missing_regions = 0
for eta_idx in range(0, nbinseta):
    for pt_idx in range(0, nbinspt):
        string_tbc = str(eta_idx) + "_" + str(pt_idx)
        local_counter = 0
        for line in enumerate(lines_good):
            if string_tbc in line[1]:
                local_counter += 1
        if local_counter < 2:
            missing_regions += 1
            print(string_tbc, "\n")
print("out of ", nbinseta*nbinspt, " regions, ", missing_regions, "are not constrained")
                
                
      

    
    
