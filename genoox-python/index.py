import requests
import gzip
from io import BytesIO
from pathlib import Path
import json
import traceback

def download_stream_zip(url):
    response = requests.get(url, stream=True,verify=False)

    if response.status_code == 200:
        with gzip.GzipFile(fileobj=BytesIO(response.content)) as gz:
            header_lines, column_line = [], None
            for line_number, line in enumerate(gz):
                line = line.decode()
                if line.startswith('##'):
                    header_lines.append(line.strip())
                elif line.startswith('#'):
                    column_line = line.strip().split('\t')
                else:
                    yield line_number,line, header_lines, column_line

def is_cache_file_exist(configurations):
    try:
        cache_file_path = Path(configurations["cache_file"])
        if cache_file_path.is_file():
            with open(cache_file_path,"r") as cach_file:
                data = json.load(cach_file)
                return data
        return {}
    except:
        print("didnt exsist cache so I creacted a file")
        
       

#GENE api request
def get_gene_from_api(chrom,pos,ref,alt ,configurations): #change cache to use file instead of memory
    payload = {
        "chr":chrom,
        "pos":pos,
        "ref":ref,
        "alt":alt,
        "reference_version":"hg19",
    }
    cache = is_cache_file_exist(configurations) #check if there is a cache file or creates it 
    
    cache_key = json.dumps(payload) #convert payload to string so i can use it as key
    
    if cache_key in cache: #check if the Gene that needed is already provided from the api
        print("Data is found in cache")
        return cache[cache_key]
    else:
    #send api request to find the gene
        try:
            response = requests.post(configurations["API_URL"],json=payload)
            if response.status_code == 200:
                response_json = response.json()
                gene_info = response_json.get("gene", "")
                cache[cache_key] = gene_info

                #save updated data to cach file
                with open(configurations["cache_file"], "w") as cach_file:
                    json.dump(cache, cach_file)
                return gene_info
            
            else:
                print(f"POST request failed with status code: {response.status_code}")

        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

#handle header_lines
def read_header_and_column_line_from_file(file):
    header_lines = []
    column_line =""
    for line in file:
        if line.startswith("##"):
            header_lines.append(line.strip())
        else:
            column_line = line.strip().split('\t')
            break
    return header_lines,column_line


def get_filtered_variant(columns, configurations): #get gene from api/cache and add to INFO column
    gene_info = get_gene_from_api(columns[configurations['chrom']], int(columns[configurations['pos']]), columns[configurations['ref']], columns[configurations['alt']], configurations) #save to variables
    columns[configurations['info']] = f"{columns[configurations['info']]};GENE={gene_info}"
    return "\t".join(columns)+"\n"

def is_relevant_sample_data(sample_data):
    return all(x != "." for x in sample_data.values())

def is_condition_apply(start,end,pos,minDP,dp):
    return start <= pos and pos <= end and dp > minDP

def handle_sample_file(sample ,samples, header_lines, samples_counter, filtered_variant, column_line):
    if sample not in samples: #new sample file is needed
        output_file_name = f"{sample}_filtered.vcf" #get smaple name and create the first sample file
        sample_column_line = column_line[:9].copy()
        sample_column_line.append(sample)
        with open(output_file_name, 'a+') as sample_file: #its add line to the end of the file 
            sample_file.writelines([line + "\n" for line in header_lines])
            sample_file.writelines("\t".join(str(x) for x in sample_column_line)+"\n")
            sample_file.writelines(filtered_variant)
            samples[sample] = output_file_name
            samples_counter[sample] =1
    else: #sample exists 
        with open(samples[sample], 'a') as sample_file:
            sample_file.writelines(filtered_variant)
            samples_counter[sample] +=1                

def save_progress(line_number,samples,samples_counter,configurations):
    progress = {}
    progress["line_number"] = line_number
    progress["samples"] = samples
    progress["samples_counter"] = samples_counter
    with open(configurations["progress_file"],"w") as progress_file:
        json.dump(progress,progress_file)
        

#handle vcf file
def process_vcf_file(start, end, minDP, limit, configurations,progress):
    if progress:
        samples = progress["samples"]
        samples_counter = progress["samples_counter"]
        progress_line = progress["line_number"]
    else:
        samples = {}
        samples_counter = {}
        progress_line = None
    try:
            for line_number,line,header_lines,column_line in download_stream_zip(configurations["data_url"]):
                sample_names = column_line[9:]  #get the samples name after the FORMAT
                if progress_line and line_number <= progress_line: #skip lines before where we left off
                    continue
                for index,sample in enumerate(sample_names): #index:sample -> 0:father, 1:mother ...
                    if sample in samples_counter and samples_counter[sample] >= limit:
                        break
                    columns = line.strip().split("\t")
                    pos = int(columns[configurations["pos"]])
                    sample_data = columns[9+index]
                    sample_data_parsed = sample_data.split(":")
                    format_data = columns[configurations['format']]
                    format_data_parsed = format_data.split(":") #split the FORMAT data column
                    sample_dict = {}
                    for i, key in enumerate(format_data_parsed): #create a dict acoorfing to FORMAT column
                        sample_dict[key] = sample_data_parsed[i]
                    if is_relevant_sample_data(sample_dict): #fillter all '.' in sample
                        if is_condition_apply(start,end,pos,minDP,int(sample_dict["DP"])): 
                            variant_with_specific_sample_only = columns[:9].copy() 
                            variant_with_specific_sample_only.append(columns[9+index]) #only spacific sample lines
                            filtered_variant = get_filtered_variant(variant_with_specific_sample_only, configurations)
                            handle_sample_file(sample ,samples, header_lines, samples_counter , filtered_variant, column_line) #if sample file exsist add relavante sample line else create
                    save_progress(line_number,samples,samples_counter,configurations)
            
    except FileNotFoundError:
        print(f"File {configurations['data_url']} is not found")
    except Exception as e:
        print("Exception occurred:")
        print("Message:", e)
        print("Stack Trace:")
        traceback.print_exc()

def read_configurations():
    try:
        with open("configurations.json","r") as file:
            configurations = json.load(file)
            if configurations:
                return configurations
            raise "Failed to load configurations!"
    except Exception as e:
        print("Exception occurred:")
        print("Message:", e)
        print("Stack Trace:")
        traceback.print_exc()

def check_saved_progress(configurations):
    try:
        print("Searching for saved progress before starting... ")
        progress_file_path = Path(configurations["progress_file"])
        if progress_file_path.is_file():
            print("loading progress file")
            with open(progress_file_path,"r") as progress_file:
                data = json.load(progress_file)
                return data
        
    except:
        print("Could not find progress file")


def main():
    configurations = read_configurations()
    progress = check_saved_progress(configurations)
    limit = configurations["limit"]
    start = configurations["start"]
    end = configurations["end"]
    minDP = configurations["minDP"]
    max_limit = configurations["max_limit"]
    if not(isinstance(limit, int) and limit < max_limit):
        print("Error: the limit must be int less than ",max_limit)
        return
    process_vcf_file(start, end, minDP, limit,configurations,progress)
    

if __name__=="__main__":
    main()