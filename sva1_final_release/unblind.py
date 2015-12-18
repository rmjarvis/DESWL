#!/usr/bin/env python
import sys
import hashlib

def get_factor(phrase="DES is blinded"):
    """
    Default code phrase for SV:

        "DES is blinded"
    
    Use the code like this:
    
        e1 = e1 / get_factor()
        e2 = e2 / get_factor()
    
    """

    #hex number derived from code phrase
    m = hashlib.md5(phrase).hexdigest()
    #convert to decimal
    s = int(m, 16)
    # last 8 digits
    f = s%100000000
    # turn 8 digit number into value between 0 and 1
    g = f*1e-8
    #get value between 0.9 and 1
    return 0.9 + 0.1*g

