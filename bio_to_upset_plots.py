from upsetplot import from_memberships
from upsetplot import plot
from matplotlib import pyplot
import re

def file2_triple(path):
  """ Prepares BIO data"""
  with open(path) as f:
    l = f.readlines()
  res= [re.split(" ", re.sub("\n", "", x)) for x in l if len(x)>2 and x[0]!="#"]
  return [x for x in res if len(x)==3]
import itertools

def get_combination(names):
  """ Returns all combination of names"""
  liste_possible = []
  for L in range(0, len(names)+1):
    for subset in itertools.combinations(names, L):
      liste_possible.append(tuple(sorted(list(subset))))
  return liste_possible

def display_res(all_res)
  for res, dic in all_res.items():
    print(res)
    for tup, nb in dic.items():
      if nb>0:
        print("  ",tup, nb)

def evaluate_predictions(data, names, i, debug):
  """ evaluates the predictions for th i{th} line"""
  truth = "P" #by default
  if data[0][i][2]=="O":
    truth = "N" #else it is a positive
  pred = {}#store results
  for j in range(len(names)):
    if data[j][i][1]==data[j][i][2]:
      this_pred = "T"
      this_truth = truth
    else:
      this_pred = "F"
      this_truth = "P" if truth =="N" else "N"#switch values
    type_res = f"{this_pred}{this_truth}"
    pred.setdefault(type_res, [])
    pred[type_res].append(names[j])
    if debug == True and type_res=="TP":
      print("-->",names[j], f"{type_res} : %s"%data[j][i])
  return pred

def get_all_res_tokens(data, names)
  liste_possible = get_combination(names)
  print(  "%i possible combinations"%len(liste_possible))
  all_res= {}
  for i in range(len(data[0])):
    pred = evaluate_predictions(data, names, i, debug)
    for type_res, liste in pred.items():
      #only make sense to get numbers from empty lists if we are on the TP case
      if len(liste)>0 or type_res=="TP":
        cle = tuple(sorted(liste))
        all_res.setdefault(type_res, {x: 0 for x in liste_possible})
        all_res[type_res][cle]+=1
  return all_res

def bio2upsetData(data, names, debug = False):
  """Transforms BIO data in upset_plot format"""
  all_res_tokens = get_all_res_tokens(data, names)
  if debug==True:
    display_res(all_res_tokens)
  return all_res_tokens

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

def verify_data(liste_files, data_in):
  """Verify the number of BIO lines and that the tokens are the same"""
  for i in range(1, len(data_in)):
    f0 = liste_files[0]
    fi = liste_files[i]
    l0 = len(data_in[0])
    li = len(data_in[i])
    assert l0==li, f"Inconsistent number of BIO lines :\n {f0} ({l0}) VS {fi} ({li})"
    for j in range(len(data_in[0])):
      tok0 = data_in[0][j][0]
      toki = data_in[i][j][0]
      assert tok0==toki, f"Not the same token line {j} : \n f{0} --> {tok0}\n f{1} --> {toki}" 

def files_2_cat(path):
  liste_files = glob.glob(f"{path}/*.tx*")
  data_in = [file2_triple(path_file) for path_file in liste_files]
  verify_data(liste_files, data_in)

  names = [re.split("/|\.", x)[-2] for x in liste_files]
  print(f"\nfilenames : {names}")
  res = bio2upsetData(data_in, names, debug=DEBUG)
  plot_graph(res, path)


#for an external usage just call the files_2_cat function with a path containg BIO files


if __name__=="__main__":
  import glob, sys
  print("Usage : python bio_to_upset_plots.py DIR")
  print("DIR is the directory where the aligned result files (BIO format) are stored")
  print("The names displayed on the graphs will be derived from the filenames")
  DEBUG = False
  if len(sys.argv)==1:
    print("\nspecify the path of the directory containg the BIO result files")
    print("\nfor example python bio_to_upset_plots.py example_data")
    print("...exiting ...")
    exit()
  elif len(sys.argv)==3:
    print("--- debug mode ---")
    DEBUG = True
  TOKENS=True
  path = sys.argv[1]
  files_2_cat(path)
