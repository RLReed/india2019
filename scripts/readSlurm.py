taskID = 2419202
Tasks = range(1, 390 + 1)

for task in Tasks:
    with open('z_trunc44-{}_{}.out'.format(taskID, task), 'r') as f:
        data = f.readlines()
        line = data[-1][:-1]
        if 'NaN' in line:
            print task, data[10][:-1], line
