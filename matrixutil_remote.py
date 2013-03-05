# -*- coding: utf-8 -*-
"""
The MIT License
===============
    
Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
'Software'), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    
"""

import collections
import scipy.sparse
from scipy import zeros,arange,mat
from scipy.sparse import vstack,hstack,csr_matrix,lil_matrix,triu
from scipy.spatial import Delaunay
from scipy.linalg import *
from pyplasm import *

try: import simplejson as json
except ImportError: import json
import requests
from numpy import *
#from StringIO import StringIO
from termcolor import colored

def csrToJSON(CSRm):
    
    ROWCOUNT = CSRm.shape[0];
    COLCOUNT = CSRm.shape[1];
    ROW = CSRm.indptr.tolist();
    COL = CSRm.indices.tolist();
    DATA = CSRm.data.tolist();
    
#   print colored('DENSE:', 'red');
#   print(CSRm.todense());
    
    JSONm = json.dumps({"ROWCOUNT":ROWCOUNT, "COLCOUNT":COLCOUNT, "ROW":ROW, "COL":COL, "DATA":DATA })

#   print colored('JSON:', 'red');
#   print(JSONm);
  
    return JSONm

def jsonToCSR(JSONm):
    
#   print colored('JSON:', 'red');
#   print(JSONm);
    
    ROWCOUNT = JSONm['ROWCOUNT'];
    COLCOUNT = JSONm['COLCOUNT'];
    ROW = JSONm['ROW'];
    COL = JSONm['COL'];
    DATA = JSONm['DATA'];

    CSRm = csr_matrix((array(DATA),array(COL),array(ROW)),shape=(ROWCOUNT,COLCOUNT));

#   print colored('DENSE:', 'red');
#   print(CSRm.todense());

    return CSRm;

def csrTranspose(CSRm):
    CSRm = CSRm.T;
    return CSRm.tocsr();

def csrProduct(CSRm1,CSRm2):
    
    JSONm1 = csrToJSON(CSRm1);
    JSONm2 = csrToJSON(CSRm2); 
    
    url = "http://cvd01.dia.uniroma3.it:3000/multiply";
    payload = {"matrixa": JSONm1, "matrixb": JSONm2};   
    req = requests.post(url, data=payload);
        
    print colored('Matrix A and matrix B sent', 'red');
    print colored('Status code: ', 'blue'), req.status_code;
    print colored('Url: ', 'blue'), req.url;
    print colored('Headers: ', 'blue'), req.headers['Content-Type'];
#   print req.content;
	
    JSONm1m2 = json.loads(req.content);
    
    print colored('Matrix A*B received', 'green');
    CSRm1m2 = jsonToCSR(JSONm1m2);
    
    return CSRm1m2
