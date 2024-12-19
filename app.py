#Student: Pierre Wensel
#Course: Bionformatics Applications-Python
#Date: January 14, 2024
#Instructor and Author of Base Application Files: Dr. A Cordomi

#This application enables a user to determine motifs, transform nucleic acid sequences, and compute properties of DNA sequences

from flask import Flask, render_template, request
from Bio.Seq import Seq
import seqtools
from Bio.SeqUtils import GC
from Bio.SeqUtils import molecular_weight

#################################################################
#
# This is my fantastic Bioinformatics Web Application
#
#################################################################


# My app
app = Flask(__name__)

@app.route("/")
def main():
    return render_template("index.html")


@app.route("/motif_input.html")
def motif_input():
    return render_template("motif_input.html")


@app.route("/motif_results.html", methods=["POST"])
def motif_results():
    input_seq = seqtools.clean(request.form["input_seq"])
    motif_seq = seqtools.clean(request.form["motif_seq"])
    input_text = seqtools.chunks(input_seq, 120)
    motif_text = seqtools.chunks(motif_seq, 120)
    if motif_seq in input_seq:
        out_text = "The sequence DOES contain the motif"
    else:
        out_text = "The sequence does NOT contain the motif"
    return render_template("motif_results.html", **locals())


@app.route("/transform_input.html")
def transform_input():
    return render_template("transform_input.html")


@app.route("/transform_results.html", methods=["POST"])
def transform_results():
    input_seq = seqtools.clean(request.form["input_seq"])
    transform = request.form["transform"]
    input_seq = Seq(input_seq)
    if transform == "translate":
        output_seq = seqtools.cutseq(input_seq).translate()
    elif transform == "transcribe":
        output_seq = input_seq.transcribe()
    elif transform == "complement":
        output_seq = input_seq.complement()
    elif transform == "reverse_complement":
        output_seq = input_seq.reverse_complement()
    output_text = seqtools.chunks(output_seq, 120)
    input_text = seqtools.chunks(input_seq, 120)
    return render_template("transform_results.html", **locals())

@app.route("/compute_input.html")
def compute_input():
    return render_template("compute_input.html")

@app.route("/compute_results.html", methods=["POST"])
def compute_results():
    input_seq = seqtools.clean(request.form["input_seq"])
    compute = request.form["compute"]
    input_seq = Seq(input_seq)
    if compute == "gc_content":
        output_seq = GC(input_seq)
    elif compute == "dna_molecular_weight":
        output_seq = molecular_weight(input_seq, "DNA")
    output_text = str(output_seq)
    #output_text = seqtools.chunks(output_seq, 120)
    input_text = seqtools.chunks(input_seq, 120)
    return render_template("compute_results.html", **locals())

# Start the app (let's keep debug=True when building the )
app.run(debug=True)

