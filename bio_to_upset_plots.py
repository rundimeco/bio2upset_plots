from upsetplot import from_memberships
from upsetplot import plot
from matplotlib import pyplot
import re
import glob

def file2_triple(path):
  """ Prepares BIO data"""
  with open(path, encoding="utf-8") as f:
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
  return sorted(liste_possible)

def display_res(all_res):
  for res, dic in all_res.items():
    print(res)
    for tup, nb in dic.items():
      if nb>0:
        print("  ",tup, nb)

def evaluate_predictions(data, names, i, debug):
  """ evaluates the predictions for th i{th} line"""
  if data[0][i][2]=="O":
    truth = "N" #else it is a positive
    pred = {"FP":[], "TN":[]}#store possible results
  else:
    truth = "P" #by default
    pred = {"TP":[], "FN":[]}#store possible results
  for j in range(len(names)):
    if data[j][i][1]==data[j][i][2]:
      this_pred = "T"
      this_truth = truth
    else:
      this_pred = "F"
      this_truth = "P" if truth =="N" else "N"#switch values
    type_res = f"{this_pred}{this_truth}"
    pred[type_res].append(names[j])
    if debug == True and type_res=="TP":
      print("-->",names[j], f"{type_res} : %s"%data[j][i])
  
  return pred

def get_all_res_tokens(data, names, debug):
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

def bio2upsetData(data, names, debug = False, verbose = False):
  """Transforms BIO data in upset_plot format"""
  all_res_tokens = get_all_res_tokens(data, names, debug)
  if verbose==True:
    display_res(all_res_tokens)
  return all_res_tokens

def write_json_file(path, content):
  """ writes some content in json format"""
  import json
  w = open(path, "w", encoding="utf-8")
  w.write(json.dumps(content, indent = 2))
  w.close()

def plot_graph(res, path, img_format = "png"):
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
    path_img = f"{path_figures}/{typ_res}.{img_format}"
    if VERBOSE==True:
      print(path_img)
    pyplot.savefig(path_img)

  print(f"--> figures stored in '{path_figures}/'")

  path_upset = f"{path}/data_upset.json"
  write_json_file(path_upset, [liste_cats, data_out])
  print(f"--> output file in upset plot format stored in '{path_upset}'")

def verify_data(list_files, data_in):
  """Verify the number of BIO lines and that the tokens are the same"""
  for i in range(1, len(data_in)):
    f0 = list_files[0]
    fi = list_files[i]
    l0 = len(data_in[0])
    li = len(data_in[i])
    assert l0==li, f"Inconsistent number of BIO lines :\n {f0} ({l0}) VS {fi} ({li})"
    for j in range(len(data_in[0])):
      tok0 = data_in[0][j][0]
      toki = data_in[i][j][0]
      assert tok0==toki, f"Not the same token line {j} : \n f{0} --> {tok0}\n f{1} --> {toki}" 

def get_names_plot(list_files, dic_names):
  """
  Gives a name for each path present in list_files
  By default: gives filename
  If dic_names contains a "translation" for each path or each filename, it gives the translation
  """
  names = [re.split("/|\.", x)[-2] for x in list_files]
  if len(dic_names)>0:
    missing_paths     = set(list_files).difference(set(dic_names.keys()))
    missing_filenames = set(names).difference(set(dic_names.keys()))
    if len(missing_paths)==0:# Translating path
      return [dic_names[x] for x in list_files]
    elif len(missing_filenames)==0:# Tranlasting filenames
      return [dic_names[x] for x in names]
  # If we are there, translation was not possible
    print("-"*20)
    print("'dic_names' is incomplete")
    print(f"It provided translations for {list(dic_names.keys())} but :")
    print(f"  - I'm missing filepath translation for {missing_paths}")
    print(f"  - I'm missing filename translation for {missing_filenames}")
    print("-"*20)
  return names

def files_2_cat(path, img_format = "png", dic_names= {}):
  list_files = glob.glob(f"{path}/*.tx*")
  names = sorted(get_names_plot(list_files, dic_names))
  print(f"\nfilenames : {names}")
  
  data_in = [file2_triple(path_file) for path_file in list_files]
  verify_data(list_files, data_in)
  res = bio2upsetData(data_in, names, debug=DEBUG, verbose = VERBOSE)
  plot_graph(res, path, img_format = img_format)

# for an external usage just call the files_2_cat function with a path containing BIO files

#TODO: vérifier DEBUG
#TODO: order pour les figures
if __name__=="__main__":
  import glob, sys
  print("Usage : python bio_to_upset_plots.py DIR")
  print("DIR is the directory where the aligned result files (BIO format) are stored")
  print("The names displayed on the graphs will be derived from the filenames")
  DEBUG = False
  VERBOSE = True
  if len(sys.argv)==1:
    print("\nspecify the path of the directory containg the BIO result files")
    print("\nfor example python bio_to_upset_plots.py example_data")
    print("... exiting ...")
    exit()
  elif len(sys.argv)==3:
    print("--- debug mode ---")
    DEBUG = True
  TOKENS=True
  path = sys.argv[1]

  ### Simple usage
  print("Simple usage : ")
  print("files_2_cat(path)")
  files_2_cat(path)

  VERBOSE = False
  ### Changing image format
  print("\n") 
  print("*"*20)
  dd= input("One can change the image format with img_format :")
  print("*"*20)
  print("files_2_cat(path, img_format=\"svg\")")
  files_2_cat(path, img_format="svg")

  ### Giving alternative names on the graphs (default: filenames)
  print("\n") 
  print("*"*20)
  dd = input("One can also give alternative names for each name in path")
  print("*"*20)
  
  # Creating dummy values for the example
  dummy_names= ["toto", "titi", "tutu"]
  liste_path = glob.glob(f"{path}/*.tx*")
  dic_path = {liste_path[i]:dummy_names[i] for i in range(len(dummy_names))}
  print(f"dic_path={dic_path}")
  print("files_2_cat(path, dic_names = dic_path)")
  files_2_cat(path, dic_names = dic_path)

  print("\n") 
  print("*"*20)
  print("What happens if dic_path is incomplete?")
  print("*"*20)

  print(f"liste path = {liste_path}")
  print(f"dic_path.keys() = {dic_path.keys()}")

  dic_path = {liste_path[i]:dummy_names[i] for i in range(2)}

  files_2_cat(path, dic_names = dic_path)
