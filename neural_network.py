from sklearn.neural_network import MLPClassifier
import csv
import argparse

def clean_data(fname):
    with open(fname, 'rb') as f:
        reader = csv.reader(f)
        out = reader.next()
        out2 = reader.next()
        outname = "row_names_{}".format(fname.split("/")[-1])
        with open(outname, 'wb') as fo:
            writer = csv.writer(fo)
            for col in out:
                writer.writerow([col])
        outname2 = "row_names_{}_2".format(fname.split("/")[-1])
        with open(outname2, 'wb') as fo:
            writer = csv.writer(fo)
            for col in out2:
                writer.writerow([col])

def split_train(data):
    return None

def split_classify(data):
    return None

def train(train_data, train_class):
    return None

def predict(model, verify_data):
    return None

def compare(prediction, verify_class):
    return None



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help="File containing pitch data", type=str)
    args = parser.parse_args()
    print(args.file)
    cleaned_data = clean_data(args.file)
    train_set, verify_set = split_train(cleaned_data)
    train_data, train_class = split_classify(train_set)
    verify_data, verify_class = split_classify(verify_set)
    model = train(train_data, train_class)
    prediction = predict(model, verify_data)
    performance = compare(prediction, verify_class)
    print performance

    
if __name__=="__main__":
    main()