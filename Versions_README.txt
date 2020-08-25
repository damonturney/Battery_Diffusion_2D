v1.1
v1.2
v1.3
v1.4









######### ALPHABETICAL LIST OF VIM COMMANDS











########## ALPHABETICAL GIT COMMANDS (MASTER LIST IS IN BACKEDUP/Software/DamonWrittenSoftware/)
cat .git/HEAD                        see what the HEAD is current pointing to
git add .                            add all changed files to the staging area
git branch                           shows the names of all local branches
git branch --all                     shows the names of all local AND remote branches
git branch --help                    show help for git branch command
git branch -v                        shows the last commit on each branch
git branch --verbose                 same as above
git branch name                      create a new branch named "name" at the HEAD of the current branch
git branch name startpoint           create a new branch named "name" at startpoint of the current branch, startpoint can be a branch name, a commit SHA (id number), or a tag name, or other ...
git checkout branchname              use "git switch branchname" instead, switches the folder contents and HEAD to be the last commit of branchname
git commit -m 'better phi function'  commit all staged files to the current branch
git init                             initiate a git repository in the current directory
git pull                             first does a fetch from the origin then integrates commits
git pull origin master               first does a fetch from the origin then integrates commits
git pull --tags
git push                             pushes all commits from the current branch to the default origin
git push --tags                      push all lightweight tags to the origin
git push origin HEAD --tags          push all commits from the current branch to origin, including lightweight tags
git push origin master               pushes committs from master to origin
git remote -v                        list the URLS of remote copies of the repository that send (fetch command) or recieve (push command) updates
git remote add                       git remote add origin https://github.com/damonturney/20190925_AlCl4_EMImCl_Battery_Simulation.git
git reset HEAD                       erases unstaged changes
git stash                            stashes unstaged changes
git status                           reports if there are uncommitted changes and if the local repository is ahead of the origin
git switch branchname                switches the folder contents and HEAD to be the last commit of branchname
git switch --help                    show help for git switch command
git tag                              show all tags on the current branch
git tag tagname                      add a lightweight tag to the current HEAD
git tag tagname commit               add a lightweight tag to a specific commit, where you identify the commit by a commit SHA (id number)
git tag -d tagname                   delete a tag locally (you need to delete it on the remote too)
git push --delete origin tagname     delete a tag on origin


SEE COMPLETE HISTORY OF THE GIT REPOSITORY
git log --all --graph --decorate --oneline --date=short --pretty



WITHIN THE .gitignore FILE
*                                    * means "everyhting", thus this line means "ignore everything"
!.gitattributes                      ! means "don't",  thus this line means "don't ignore .gitattributes"
!.gitignore                          means "don't ignore .gitignore"
!*.py                                don't ignore all files that end with .py
!versions.txt                        don't ignore versions.txt
!*/                                  don't ignore all subdirectories, but inside each subdirectory all other gitignore rules are applied

 


SETUP GIT ON LOCAL COMPUTER AND ENABLE TRANSFERS TO GITHUB WITHOUT PASSWORDS
git config --global user.name "Damon Turney"
git config --global user.email "damonturney@gmail.com"
ssh-keygen -f ~/.ssh/id_rsa -t rsa -C "damonturney@gmail.com"   And enter a password when it prompts
copy the rsa key from ~/.ssh/id_rsa.pub with pbcopy < ~/.ssh/id_rsa.pub and paste it into the SSH key entry in GitHub webaccount "Settings --> SSH and GPG keys"



INITIATE A GIT REPOSITORY ON LOCAL COMPUTER
git init
vim .gitignore
git add .
git commit -m "first commit"
git tag v0.0
git remote add origin https://github.com/damonturney/20190925_AlCl4_EMImCl_Battery_Simulation.git



COMMIT CHANGES THEN PUSH CHANGES TO REMOTE
git add .                             add all changed files to the staging area
git commit -m 'better phi function'   commit all staged files to the current branch
git push origin HEAD --tags           push all commits from current branch to origin, including lightweight tags



VERSIONING
git tag tagname commit               add a lightweight tag to a specific commit, where you identify the commit by a commit SHA (id number)


PUSH ALL BRANCHES AND ALL TAGS TO REMOTE
git push --all
git push --tags


CREATE A BRANCH THEN SWITCH TO IT
git branch name                        create a new branch named "name" at the HEAD of the current branch
git switch name                        switches the folder contents and HEAD to be the last commit of branchname
