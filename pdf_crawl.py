import os
import csv
import json

#从导出的csv文件中提取出物相的pdf卡片号，分子式，晶胞参数，并存入json文件中。
json_data = []
with open(r"D:\03-CS\ED-TOOL\json_data.json", "w") as jd:
    path = r"D:\03-CS\ED-TOOL\cell parameter\\"
    file_list = os.listdir(path)
    print(file_list)
    for file in file_list:
        with open(path+file) as f:
            f_csv = csv.reader(f)
            headers = next(f_csv)
            for row in f_csv:
                #row[2]: 01-234-5678, we just need 34-5678
                card_number = row[2][4:]
                name = row[3]
                molecular_formula = ''.join(row[4].split(' '))
                raw_cell_parameter = row[7].split('-')
                cell_parameter = {}
                for i in raw_cell_parameter:
                    j = i.strip().split(' ')
                    cell_parameter[j[0]] = float(j[1])
                a = cell_parameter['a']
                b = cell_parameter['b']
                c = cell_parameter['c']
                #有些结构没有直接给出角度的值
                alpha = cell_parameter.get('alpha', 90.00)
                beta = cell_parameter.get('beta', 90.00)
                if 'Hexagonal' in file:
                    gamma = cell_parameter.get('gamma', 120.00)
                else:
                    gamma = cell_parameter.get('gamma', 90.00)
                data = dict(card_number = card_number, compound_name = name, molecular_formula = molecular_formula,
                            a = a, b = b, c = c, alpha = alpha, beta = beta, gamma = gamma)
                print(gamma)
                json_data.append(data)
    result = json.dump(json_data, jd, indent=2)   #加indent，输出格式化的json， 有换行和缩进，好看！

print('finished')

