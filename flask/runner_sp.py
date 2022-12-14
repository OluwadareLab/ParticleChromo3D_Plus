def make_tree(path):
    tree = dict(name=os.path.basename(path), children=[])
    try: lst = os.listdir(path)
    except OSError:
        pass #ignore errors
    else:
        for name in lst:
            fn = os.path.join(path, name)
            if os.path.isdir(fn):
                tree['children'].append(make_tree(fn))
            else:
                tree['children'].append(dict(name=name))
    return tree

from flask import Flask, render_template, jsonify , request, send_file
import subprocess
from werkzeug.utils import secure_filename
import os

import glob

app = Flask(__name__)

@app.route("/process")
def home():
  # args
  ifname = request.args.get('ifname', type=str)
  ss = request.args.get('ss', type=str)
  itt = request.args.get('itt', type=str)
  threshold = request.args.get('threshold', type=str)
  randRange = request.args.get('randRange', type=str)
  lf = request.args.get('lf', type=str)
  outFile = request.args.get('outFile', type=str)
  email = request.args.get('email', type=str)
  # Custom upload
  if (ifname == "uploaded_file"):
    list_of_files = glob.glob('/apt/repo/*')
    ifname = os.path.basename(max(list_of_files, key=os.path.getctime))

 
  # cmd
  thecmd = "python3 Ps.py -ss " + ss \
    +  " -itt " + itt \
    +  " --threshold " + threshold \
    +  " --randRange " + randRange \
    +  " -lf " + lf \
    +  " --outfile " + "out/"+outFile \
    +  " --email " + email \
    + " /apt/upload/" + ifname
  if (len(email) < 2):
    return "ERROR: Email required"
  else:
    p = subprocess.Popen(thecmd, shell=True, stdout=subprocess.PIPE, 
                     stderr=subprocess.PIPE)
    return "Processing has begun. You will recieve an email upon its completion."  
    #return subprocess.check_output(thecmd, shell=True) 
  #return subprocess.check_output("python3 Ps.py chr23_matrix.txt", shell=True) 

@app.route("/neo")
def neo():
  return subprocess.check_output("python3 pdbToNeo4j.py", shell=True)



@app.route('/form')
def form():
    return render_template('form.html')
 
@app.route('/upload', methods = ['POST', 'GET'])
def upload():
    if request.method == 'POST':
        f = request.files['file']
        f.save("//apt//repo//"+secure_filename(f.filename))
        return "File saved successfully"

@app.route('/uploaded')
def dirtree():
    path = '//apt//repo//'
    return render_template('dirtree.html', tree=make_tree(path))

@app.route('/download')
def download_file():
    ofname = request.args.get('ofname', type=str)

    path = "/apt/" + ofname + ".pdb"
    return send_file(path, as_attachment=True)

if __name__ == "__main__":
    app.run(host='0.0.0.0',port = 5001,debug=True)
