#!/bin/sh

rm -rf docs/api
epydoc --html -v -o api gx

