# learngit
For learning use [Reference link](https://www.liaoxuefeng.com/wiki/896043488029600 "廖雪峰的官方网站-Git教程")
## Command line tools
### 1. Check git verion:
```
$ git version
git version 2.30.1 (Apple Git-130)
```      
### 2. Building local repository
Creat a new directory
```
$ mkdir learngit
$ cd learngit
$ pwd
/Users/zhuyuehua/learngit
```
Use `git init` to establish a Git repository
```
$ git init
Initialized empty Git repository in /Users/michael/learngit/.git/
```
### 3. Add file to repository
```
$ git add readme.txt
$ git commit -m "wrote a readme file"
```
Commit several files at once
```
$ git add file1.txt
$ git add file2.txt file3.txt
$ git commit -m "add 3 files."
```
### 4. Remove file from repository
Remove from local
```
$ rm test.txt
```
Remove from Git and commit
```
$ git rm test.txt
rm 'test.txt'

$ git commit -m "remove test.txt"
[master d46f35e] remove test.txt
 1 file changed, 1 deletion(-)
 delete mode 100644 test.txt
```
Restore local file from repository (版本库里的版本替换工作区的版本)
```
$ git checkout -- test.txt
```
