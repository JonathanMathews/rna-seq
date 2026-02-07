#!/bin/bash
#get_fastq_sra.sh

usage() {
  echo "-h Help documentation for get_from_refseq.sh"
  echo "-i  SRA Accession Run ID (SRR or ERR), ProjectID (PRJ) or SampleID (SAM)"
  echo "Example: bash getSraGeo.sh -i ${id}"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :i:f:h opt
do
    case $opt in
        i) id=$OPTARG;;
        f) format=$OPTARG;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }

shift $(($OPTIND -1))

if [[ -z $id ]]
then
    usage
fi

if [[ -f ${id}.design.txt ]]
then
    rm ${id}.design.txt
fi

isSRA=$(echo $id | grep 'CRS\|CRX\|DRP\|DRR\|DRS\|DRX\|ERP\|ERR\|ERS\|ERX\|SRP\|SRR\|SRS\|SRX\|SAMD\|SAMEA\|SAMEG\|SAMN\|PRJ')

writehtml () {
    echo "<html>
 <head>
   <title>
   Sequence Not Found
   </title>
 </head>
 
 <body>
   ${1}
 </body>
 </html>" > err.html
}
fields="run_accession,sample_accession,sample_alias,sample_title,cell_type,cell_line,secondary_sample_accession,tissue_type,bio_material,experimental_factor"

if [[ $id =~ ^GS ]]
then
echo $id
    ffq -o out.json -l 1 ${id}
    urls=$(cat out.json | jq -r '.[] | .supplementary_files' | jq -r '.[] | .url')
    for j in $urls
    do
	curl -LO ${j}
    done		
elif [[ -n $isSRA ]]
then
    if [[ $id =~ ^SRR ]] || [[ $id =~ ^ERR ]]
    then
	srrids=$id
    else 
	curl -o query.txt "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${id}&result=read_run"
	srrids=$(cut -f 1 query.txt |grep -v run_accession)
    fi
    if [[ "$srrids" =~ .*"not valid for search requests on read_run data.".* ]] || [[ -z "$srrids" ]]; then
        writehtml "No sequences found for entry ${id}."
        rm *.design.txt
        exit 0
    else
        for i in $srrids
        do
        ffq -o ${i}.json $i
        fqurls=$(cat ${i}.json | jq -r '.[] | .files.ftp' | jq -r '.[] | .url' |grep 'fq\|fastq')
        #Note gcp is more complex since it needs authentication
        #gcp=$(cat ${i}.json | jq -r '.[] | .files.gcp' | jq -r '.[] | .url')
        sraurl=$(cat ${i}.json | jq -r '.[] | .files.ncbi' | jq -r '.[] | .url')
        curl -o run_info.${i}.txt -X GET "https://www.ebi.ac.uk/ena/portal/api/search?dataPortal=ena&query=run_accession=${i}&result=read_run&fields=${fields}"
        grep -v accession run_info.${i}.txt >> ${id}.design.txt
        if [[ -n $fqurls ]] && [[ $fqurls != 'null' ]]
        then
            for j in $fqurls
            do
            curl -LO ${j}
            done
        else
            for j in $sraurl
            do
            curl -LO ${j}
            sra=$(basename ${j})
            fasterq-dump $sra
            done
        fi
        done
        count=`ls -1 *fastq 2>/dev/null | wc -l`
        if [ $count != 0 ]
        then
        for j in *fastq
        do
        gzip $j
        done
        fi 
    fi
fi
