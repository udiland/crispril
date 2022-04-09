from flask import Flask, flash, request, redirect, url_for, send_from_directory
from werkzeug.utils import secure_filename
from crispys_code import Stage0
from pydantic import BaseModel
from pydantic_webargs import webargs
import os

app = Flask(__name__)


UPLOAD_FOLDER = '/input'
ALLOWED_EXTENSIONS = {'txt', 'fa', 'fasta'}
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
    """
    function to check the input file extention
    """
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route("/upload", methods=['GET', 'POST'])
def upload_file():
    """
    this function check properties of the input file and saves it in /input folder
    """
    if request.method == 'POST':
        # check if the post request has the file part
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        # If the user does not select a file, the browser submits an
        # empty file without a filename.
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

    return '''
        <!doctype html>
        <title>Upload new File</title>
        <h1>Upload new File</h1>
        <form method=post enctype=multipart/form-data>
          <input type=file name=file>
          <input type=submit value=Upload>
        </form>
        '''

    # os.system("ls -lh /input")


class RequestBase(BaseModel):
    fasta_file: str
    alg: str
    where_in_gene: float
    Omega: float
    df_targets: str
    internal_node_candidates: int
    PS_number: int


@app.route("/")
def hello():
    return "Hello World from Flask!"


@app.route('/download', methods=["GET"])
def get_output():
    """Download crispys result file to /crispys_out folder."""
    try:
        return send_from_directory("/crispys_out/", "output.csv", as_attachment=True)
    except FileNotFoundError:
        abort(404)


@app.route("/crispys", methods=["POST"])
@webargs(body=RequestBase)
def run_crispys(**kwargs):
    """
    Get crispys arguments from user
    """
    print(kwargs)
    payload = kwargs["payload"]
    crispys_req = RequestBase(
        fasta_file=payload["fasta_file"],
        alg=payload["alg"],
        where_in_gene=payload["where_in_gene"],
        Omega=payload["Omega"],
        df_targets=payload["df_targets"],
        internal_node_candidates=payload["internal_node_candidates"],
        PS_number=payload["PS_number"]
    )
    # order the parameters in a list
    crispys_args = [
        crispys_req.fasta_file,
        crispys_req.alg,
        crispys_req.where_in_gene,
        crispys_req.Omega,
        crispys_req.df_targets,
        crispys_req.internal_node_candidates,
        crispys_req.PS_number
    ]

    # print(crispys_args)


    # get the input file path
    input_fasta = "/input/" + crispys_args[0]
    # run crispys
    Stage0.CRISPys_main(input_fasta, "/crispys_out", crispys_args[1], crispys_args[2],
                        1, crispys_args[3], crispys_args[4], "outfile", 20,
                        20, 0, crispys_args[5], crispys_args[6])

    # check there is an output
    os.system( "wc -l /crispys_out/output.csv" )


    return {"args":crispys_args}

if __name__ == "__main__":
    # Only for debugging while developing
    app.run(host="0.0.0.0", debug=True, port=80)