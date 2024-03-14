with open('passed_regions.txt', 'r') as f:
    lines = f.readlines() # read and store all lines into list, note list index starts from 0

lines_good = []
with open('passed_regions_edit.txt', 'w') as fout: # print in a file lines that are good for the next step of the check
    for line in enumerate(lines): # enumerate because I need to compare the order of lines
        # for label a_b_c_d it gives c_d_a_b
        rev = line[1].split("_")[2] +"_"+ line[1].replace('\n', '').split("_")[3] +"_"+ line[1].split("_")[0] +"_"+ line[1].split("_")[1] 
        flag = 0
        for line_tbc in enumerate(lines):
            # compare each line in input file with its reverse
            if line_tbc[1].replace('\n', '') == rev: # if both the line and its reverse exist, keep only the line
                flag = 1 
                if line[0] < line_tbc[0]: # make sure you don't save the reverse as well at a later iteration
                    lines_good.append(line[1].replace('\n', ''))
                    print(line[1].replace('\n', ''), file=fout)
        if flag == 0: # if the line doesn't have a reverse, keep it
            lines_good.append(line[1].replace('\n', ''))
            print(line[1].replace('\n', ''), file=fout)

#generate the IDs for all individual muon eta_pt regions e.g 0_0, 0_1 etc. and check if there are any IDs that appear less than twice
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
            print(string_tbc, ": not constrained", "\n")
print("out of ", nbinseta*nbinspt, " regions, ", missing_regions, "are not constrained")
                
                
      

    
    
