#!/usr/bin/env python

"""
Created on Fri Jun 17 08:36:09 2016

@author: charlesgulian
"""

# Learning argparse: fibn.py
import argparse

def fib(n): # Calculate the nth Fibonacci #
    a,b = 0,1
    for i in range(n):
        a,b = b,a+b
    return a
    
def Main():
    parser = argparse.ArgumentParser() # Create argument parser
    # Then "add argument:" I think this works sort of like a dictionary, where "num"
    # is the name of the argumen and type=() specifies type of argument

    # "Add a group" (???)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-v','--verbose',action='store_true')
    group.add_argument('-q','--quiet',action='store_true')
    
    parser.add_argument('num',help='Enter the index of the Fibonacci number you wish to calculate',type=int)
    parser.add_argument('-o','--output',help='Output result to .txt file',action='store_true')    
    args = parser.parse_args()    
    
    result = fib(args.num)
    if args.verbose:
        print 'The '+str(args.num)+'th Fibonacci number is '+str(result)
    elif args.quiet:
        print str(result)
    else:
        print 'fib('+str(args.num)+') = '+str(result)
    
    if args.output:
        f = open('fibo.txt','a')
        f.write(str(result)+'\n')

if __name__ == '__main__': # What does this do?
    Main()
    
