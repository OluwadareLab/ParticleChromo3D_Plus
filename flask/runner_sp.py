
## imports
from flask import Flask, render_template, jsonify , request, send_file, redirect
import subprocess
from werkzeug.utils import secure_filename
import os

import re
import glob

# constants
regex = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,7}\b'

# Non flask functions
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

# https://www.geeksforgeeks.org/check-if-email-address-valid-or-not-in-python/
def check(email):
 
    # pass the regular expression
    # and the string into the fullmatch() method
    return re.fullmatch(regex, email)
        

## Flask
app = Flask(__name__)

@app.route("/process")
def home():
  # args
  ifname = "/apt/repo/" + request.args.get('ifname', type=str)
  ss = request.args.get('ss', type=str)
  itt = request.args.get('itt', type=str)
  threshold = request.args.get('threshold', type=str)
  randRange = request.args.get('randRange', type=str)
  lf = request.args.get('lf', type=str)
  outFile = request.args.get('outFile', type=str)
  email = request.args.get('email', type=str)
  # Custom upload
  if (ifname == "/apt/repo/uploaded_file"):
    list_of_files = glob.glob('/apt/upload/*')
    ifname = "/apt/upload/" + os.path.basename(max(list_of_files, key=os.path.getctime))

 
  # cmd
  thecmd = "python3 Ps.py -ss " + ss \
    +  " -itt " + itt \
    +  " --threshold " + threshold \
    +  " --randRange " + randRange \
    +  " -lf " + lf \
    +  " --outfile " + "out/"+outFile \
    +  " --email " + email \
    +  " " + ifname
  if (not check(email)):
    return "ERROR: Valid Email Required <br> Use back to return home"
  else:
    p = subprocess.Popen(thecmd, shell=True, stdout=subprocess.PIPE, 
                     stderr=subprocess.PIPE)
    return "Processing has begun. You will recieve an email upon its completion."  
    #return subprocess.check_output(thecmd, shell=True) # uncomment this to debug the python

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
        f.save("/apt/upload/"+secure_filename(f.filename))
        return redirect("http://" + os.environ['HOSTNAME_BE'] +":8080/ParticleChromo3D/main.html", code=302)
        return "File saved successfully"

@app.route('/uploaded')
def dirtree():
    path = '/apt/upload/'
    return render_template('dirtree.html', tree=make_tree(path))

@app.route('/download')
def download_file():
    ofname = request.args.get('ofname', type=str)

    path = "/apt/" + ofname + ".pdb"
    return send_file(path, as_attachment=True)

@app.route('/convert')
def convert():
  iffname = request.args.get('filename', type=str)
  path = "//apt//upload//" + iffname
  thecmd = "python3 TransformVCM.py " + path
  p = subprocess.Popen(thecmd, shell=True, stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)

  path2="/apt/convert.out"
  return send_file(path2, as_attachment=True)


if __name__ == "__main__":
    app.run(host='0.0.0.0',port = 5001,debug=True)
