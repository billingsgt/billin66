#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Project                 Honors
Name:                   Grant Billings
Time to completion:     15h
Comments:               This program interfaces with Ensembl Genomes' REST API
                        to pull down gene models from functional annotations
                        for further analysis.  Genes can be sorted based on
                        their broad classification, and then are accessible
                        as CSV files

Sources:                REST API documentation
                        Py Documentation for sqlite3, requests, json, and csv
"""

# used for building the database to query for chromosome placement of genes
import sqlite3
# used for querying REST API
import requests
# used for storing and interpreting the data from REST API
import json
# use for checking for files in folder
import os
# used for saving the end products of the program
import csv

# the root of the server to which all requests are sent
server = "http://rest.ensemblgenomes.org"


def welcome():
    """Prints the welcome message"""
    print("This program allows you to interface with Ensembl Plants.")
    print()
    print("It utilizes the REST API to pull content from the service.")
    print("Right now, you can dowload data from a single species at a time.")
    print()


def choose_species():
    """Gets user input to choose the species for further analysis."""
    def get_all_species():
        """GET info/species"""
        ext = "/info/species?division=EnsemblPlants"
        r = requests.get(server+ext, headers={"Content-Type":
                                              "application/json"})
        allDict = dict()
        for taxon in json.loads(r.text)["species"]:
            allDict[taxon["display_name"]] = taxon
        return allDict

    print("Here is a list of currently available species with genomes:")
    print()

    allMap = get_all_species()  # function call to the GET here

    # menu options; done multiple times
    options = {str(i): n for i, n in enumerate(sorted(allMap.keys()), 1)}

    # print the menu
    for i, n in options.items():
        print("\t{0:2} {1}".format(i, n))

    opt = input("Enter the number for the species you are interested in: ")

    while opt not in options.keys():
        opt = input("We couldn't find that one. Try again: ")

    print()

    # returns the information for the chosen species
    return allMap[options[opt]]


def choose_genome(_name):
    """
    Future implementations will allow the user to choose of one of any
    genomes.
    The procedure I was doing for this was not working since one of the
    commands stopped working (for some reason)
    """
    ext = "/info/assembly/" + _name

    r = requests.get(server+ext, headers={"Content-Type":
                                          "application/json"})

    genomeInfo = json.loads(r.text)

    print("We will use the following assembly: {0}".format(genomeInfo[
            "assembly_name"]))

    # returns all of the information from the genome assembly to be used
    return genomeInfo


def choose_gene_type(_name):
    """Allows the user to choose a specific type of gene for analysis."""
    ext = "/info/biotypes/" + _name

    r = requests.get(server+ext, headers={"Content-Type":
                                          "application/json"})

    # sort all of the biotypes available for this annotation
    annotationInfo = sorted(json.loads(r.text),
                            key=lambda d: d["biotype"].lower())

    options = {str(i): n for i, n in enumerate(sorted([d["biotype"] for
               d in annotationInfo if "gene" in d["objects"]]), 1)}

    print()
    for i, n in options.items():
        print("\t{0:2} {1}".format(i, n))

    opt = input("What type of annotations are you interested in? Select one: ")

    while opt not in options.keys():
        opt = input("We couldn't find that one. Try again: ")

    # return the data for the type of annotation that the user choses
    return [d for d in annotationInfo if d["biotype"] == options[opt]][0]


def get_gene_models(_name, _type):
    """Gets all of the gene models and saves them in JSON format."""
    print()

    allOutfile = _name + "_ALL.txt"

    if allOutfile in os.listdir():
        print("Loading from the file in the directory.")
        print()
        # to load from a previous file:
        geneModels_all = json.load(open(allOutfile, "r"))
    else:
        # query the server for the data
        print("Gene models will be saved in the file {0}".format(allOutfile))
        print()
        print("Downloading gene models now. This may take a few seconds.")
        print("All of the gene models for this taxa will be downloaded.")
        print()
        ext = "/lookup/genome/" + _name
        r = requests.get(server+ext, headers={"Content-Type":
                                              "application/json"})
        geneModels_all = json.loads(r.text)
        with open(allOutfile, "w") as outfile:
                json.dump(geneModels_all, outfile)

    print("There are a total of {0} gene models in the unfiltered list."
          .format(len(geneModels_all)))

    geneModels_filtered = [d for d in geneModels_all if d["biotype"] == _type]

    print("There is a total of {0} gene models of type {1} in filtered list."
          .format(len(geneModels_filtered), _type))
    print()

    filteredOutfile = _name + "_" + _type + ".txt"

    with open(filteredOutfile, "w") as outfile:
        json.dump(geneModels_filtered, outfile)

    print("They are saved in the file {0} in json format."
          .format(filteredOutfile))

    return geneModels_filtered


def run_DB(genes, karyotype, _name, bt):
    """Creates a database, queries it for the requested gene models."""
    print()
    print("Building a database for the gene models in memory...")
    conn = sqlite3.connect(":memory:")
    c = conn.cursor()

    headers = ", ".join([i + " text" for i in genes[0].keys()
                         if i not in ["biotype", "coord_system",
                                      "seq_region_synonyms"]])

    c.execute("CREATE TABLE gene_models ({0})".format(headers))

    headers = tuple([i for i in genes[0].keys()
                     if i not in ["biotype", "coord_system",
                                  "seq_region_synonyms"]])

    for gene in genes:
        args = tuple([gene[i] if i in gene.keys() else None
                      for i in headers])
        c.execute("INSERT INTO gene_models VALUES ({0})"
                  .format(",".join(["?"] * len(headers))), args)

    conn.commit()

    print()
    print("Database constructed in memory.")
    print()
    print("Here, you can filter the results and print them to a file.")
    print()
    print("At any time, type 'q' in order to exit.")
    print("For now, the data will be written to a CSV file.")
    print()

    options = {str(i): n for i, n in enumerate(karyotype, 1)}
    options[str(len(options) + 1)] = "other_scaffolds"

    while True:
        print("Here are the regions/chromosomes to choose from: ")
        print()

        for k, v in options.items():
            print("\t{0:>2} {1}".format(k, v))

        _chr = input("Select a number from above: ")

        while _chr not in options.keys() and _chr.lower() != "q":
            _chr = input("Not a valid choice. Try again: ")

        print()

        if _chr.lower() == "q":
            break

        _chr = options[_chr]
        _chr = tuple([_chr])

        if _chr != ("other_scaffolds",):
            c.execute("""SELECT * from gene_models WHERE seq_region_name = ?
                      """, _chr)
            data = c.fetchall()
            conn.commit()
        else:
            c.execute("""SELECT * from gene_models WHERE seq_region_name
                      NOT IN ({})"""
                      .format(",".join(["?"] * (len(options) - 1))),
                      karyotype)
            data = c.fetchall()
            conn.commit()

        if not data:
            print("No genes of this type found in this region.")
            continue

        write_to_csv(data, headers, _name, bt, *_chr)

    conn.close()


def write_to_csv(data, header, _name, bt, _chr):
    """Writes the query to a CSV file."""
    outFileName = "_".join([_name, bt, "chr", _chr]) + ".txt"

    print("Outputing as CSV into {0}...".format(outFileName))

    with open(outFileName, "w") as csvfile:
        myWriter = csv.writer(csvfile, delimiter=",")
        myWriter.writerow(header)
        myWriter.writerows(data)

    print("Data written.")
    print()


def main():
    """The main loop."""
    welcome()
    mySpeciesInfo = choose_species()
    genomeChoice = choose_genome(mySpeciesInfo["name"])
    karyotype = tuple(genomeChoice["karyotype"])

    annotationType = choose_gene_type(mySpeciesInfo["name"])

    print()
    print("Okay, so we're using the assembly <{0}> of {1} for {2}{3}."
          .format(genomeChoice["assembly_name"],
                  mySpeciesInfo["display_name"],
                  annotationType["biotype"],
                  " genes"))

    geneModels = get_gene_models(mySpeciesInfo["name"],
                                 annotationType["biotype"])

    run_DB(geneModels, karyotype,
           mySpeciesInfo["name"], annotationType["biotype"])

    print("Exiting program.")


if __name__ == "__main__":
    main()
