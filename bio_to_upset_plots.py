from upsetplot import from_memberships
from upsetplot import plot
from matplotlib import pyplot
import re

def file2_triple(path):
  """ Prepares BIO data"""
  f = open(path)
  l = f.readlines()
  f.close()
  res= [re.split(" ", re.sub("\n", "", x)) for x in l if len(x)>2]
  return [x for x in res if len(x)==3]
import itertools

def get_combination(names):
  """ Returns all combination of names"""
  liste_possible = []
  for L in range(0, len(names)+1):
    for subset in itertools.combinations(names, L):
      liste_possible.append(tuple(sorted(list(subset))))
  return liste_possible

def bio2upsetData(data, names):
  """Transforms BIO data in upset_plot format"""
  liste_possible = get_combination(names)
  print(  "%i possible combinations"%len(liste_possible))
  all_res = {}
  for i in range(len(data[0])):

    truth = "P" #by default
    if data[0][i][1]=="O":
      truth = "N" #else it is a negative

    pred = {"T":[], "F":[]}#to store True and False predictions
    for j in range(len(names)):
      if data[j][i][1]==data[j][i][2]:
        pred["T"].append(names[j])
      else:
        pred["F"].append(names[j])

    for p, liste in pred.items():
      type_res = f"{p}{truth}"
      #only make sense to get numbers from empty lists if we are on the TP case
      if len(liste)>0 or type_res=="TP":
        cle = tuple(sorted(liste))
        all_res.setdefault(type_res, {x: 0 for x in liste_possible})
        all_res[type_res][cle]+=1
  return all_res

def write_json_file(path, content):
  """ writes some content in json format"""
  import json
  w = open(path, "w")
  w.write(json.dumps(content, indent = 2))
  w.close()

def plot_graph(res, path):
  """ From upset_plot data Plot upset plots and store corresponding data"""

  path_figures = f"{path}/figures"
  import os
  os.makedirs(path_figures, exist_ok = True)

  for typ_res, dic in res.items():
    liste_cats = sorted(dic.keys())
    data_out = []
    for cat in liste_cats:
      data_out.append(dic[cat])
    example = from_memberships(liste_cats, data = data_out)
    plot(example)
    pyplot.savefig(f"{path_figures}/{typ_res}.png")

  print(f"  figures stored in '{path_figures}/'")

  path_upset = f"{path}/data_upset.json"
  write_json_file(path_upset, [liste_cats, data_out])
  print(f"  output file in upset plot format stored in '{path_upset}'")

def files_2_cat(path):
  liste_files = glob.glob(f"{path}/*.tx*")
  data_in = [file2_triple(path_file) for path_file in liste_files]
  names = [re.split("/|\.", x)[-2] for x in liste_files]
  print(f"\nfilenames : {names}")
  res = bio2upsetData(data_in, names)
  plot_graph(res, path)


#for an external usage just call the files_2_cat function with a path containg BIO files


if __name__=="__main__":
  import glob, sys
  print("Usage : python bio_to_upset_plots.py DIR")
  print("DIR is the directory where the aligned result files (BIO format) are stored")
  print("The names displayed on the graphs will be derived from the filenames")
  
  if len(sys.argv)!=2:
    print("\nspecify the path of the directory containg the BIO result files")
    print("\nfor example python bio_to_upset_plots.py example_data")
    print("...exiting ...")
    exit()
  
  path = sys.argv[1]
  files_2_cat(path)
