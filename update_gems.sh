#!/bin/bash

while read line
do
  gem install $line
done < $1
