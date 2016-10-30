
from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfVectorizer
import numpy as np

from nltk import word_tokenize          
from nltk.stem import WordNetLemmatizer 
class LemmaTokenizer(object):
     def __init__(self):
         self.wnl = WordNetLemmatizer()
     def __call__(self, doc):
         return [self.wnl.lemmatize(t) for t in word_tokenize(doc)]
         


#==============================================================================
#All the  categories
 
#['alt.atheism',
# 'comp.graphics',
# 'comp.os.ms-windows.misc',
# 'comp.sys.ibm.pc.hardware',
# 'comp.sys.mac.hardware',
# 'comp.windows.x',
# 'misc.forsale',
# 'rec.autos',
# 'rec.motorcycles',
# 'rec.sport.baseball',
# 'rec.sport.hockey',
# 'sci.crypt',
# 'sci.electronics',
# 'sci.med',
# 'sci.space',
# 'soc.religion.christian',
# 'talk.politics.guns',
# 'talk.politics.mideast',
# 'talk.politics.misc',
# 'talk.religion.misc']
#==============================================================================

 
#categories = ['comp.graphics','comp.os.ms-windows.misc','comp.sys.ibm.pc.hardware','comp.sys.mac.hardware','comp.windows.x']
# 

#categories =[  'rec.autos',
# 'rec.motorcycles',
# 'rec.sport.baseball',
# 'rec.sport.hockey']

#==============================================================================
#fetch_20newsgroups data
#categories = ['sci.crypt','sci.electronics','sci.med','sci.space']

categories =[  'rec.autos','rec.motorcycles','rec.sport.baseball','rec.sport.hockey']
newsgroup = fetch_20newsgroups(subset = 'train', categories = categories, shuffle=True, random_state = 1) 
list(newsgroup.target_names)
list(newsgroup.target)
cvectorizer = CountVectorizer(tokenizer=LemmaTokenizer(),stop_words='english')
tfvectorizer = TfidfVectorizer(tokenizer=LemmaTokenizer(),stop_words='english')


###word counts
count = cvectorizer.fit_transform(newsgroup['data']).toarray()
cvectorizer.vocabulary_

count.shape

###tf_idf vector
tf_idf = tfvectorizer.fit_transform(newsgroup['data']).toarray() 
tfvectorizer.vocabulary_


#==============================================================================
#Sort the average tf_idf in a descending order, selet top k1 words with largest tf-idf
#k1 = 2000
k1 = 2000
avg_tf_idf = tf_idf.mean(0)
words_index1 = np.array(np.argsort(-avg_tf_idf)).flatten().tolist()[0:k1]

words_index1
#avg_count = count.mean(0)
#words_index1 = np.array(np.argsort(-avg_count)).flatten().tolist()[0:k1]


#==============================================================================
#Check each selected word and neglect  the words containing non-alphabetical characters,e.g
#numbers, slash or punctuation
words_index2 = []
cdict = cvectorizer.vocabulary_
for i in words_index1:
##Be very careful!! Dict is un-ordered, we need to get the key corresponding to a given value!
    word =  list(cdict.keys())[list(cdict.values()).index(i)]
    if word.isalpha():
        words_index2.append(i)        

print (len(words_index2))


#==============================================================================
#Select k variables randomly for each of the gaussian,  count and binary category
#np.random.seed(1001)
np.random.seed(1)
#k = 30
k = 50
randomlist = np.random.choice(range(len(words_index2)),3*k,replace=False)
words_index = [words_index2[i] for i in randomlist]



words_chosen = [list(cdict.keys())[list(cdict.values()).index(i)] for i in words_index]
print (words_chosen)


###create arries
array_gaussian = tf_idf[:,words_index[0:k]]
array_count = count[:,words_index[k:2*k]]
array_binary = 1*(count[:,words_index[2*k:3*k]]>0)

data = np.column_stack((array_gaussian,array_count,array_binary))




#np.savetxt('C:/Users/student/Dropbox/iGGL/code/refs/mixed_data/newsgroup/data.txt',data, fmt='%.3f')

#file = open('C:/Users/student/Dropbox/iGGL/code/refs/mixed_data/newsgroup/varnames.txt','w')
#file.writelines(['%s\n' % item  for item in words_chosen])
#file.close()

#==============================================================================
#Show the documents containing a given  keyword 
def show_doc(keyword, data, k,m):
    index_vector = list(np.nonzero(count[:,cdict.get(keyword)])[0])[k:m]
    for i in index_vector:
        print (data[i])
        print (i)

show_doc('callison', newsgroup['data'], 10,20 )

cvectorizer.fit_transform(['US local']).toarray()
cvectorizer.vocabulary_


####Check basics
count[:,words_index[0:(3*k)]]
tf_idf[:,words_index[0:(3*k)]]
print(np.amax(sum(count[:,words_index[k:(2*k)]])))
print(np.amin(sum(count[:,words_index[k:(2*k)]])))
print(np.amax(sum(count[:,words_index1])))
print(np.amin(sum(count[:,words_index1])))
