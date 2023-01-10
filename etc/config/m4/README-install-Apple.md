# Version 4.12 of FreeFEM 
==============================================

Be carefull  the Macos security (com.apple.quarantine)

To install FreeFEM
- mount the dmg file 
- Launch a new terminal windows and  in dmg volume do: 
   remark: the terminal app is  menu Finder Go->Utility->Terminal 
 
**to install and remove com.apple.quarantine  when you copy  FreeFem++.app in /Application directory   do :**
 
>  cd /Volumes/FreeFEM-4.12-Apple-M1
>  bash ./Install-app.sh

For Autorisation see [Safely open apps on your Mac](https://support.apple.com/en-us/HT202491) or 
[Ouvrir des apps en toute sécurité sur votre Mac](https://support.apple.com/fr-fr/HT202491)

After try (cut and paste):
> /Applications/FreeFem++.app/Contents/ff-4.12/bin/FreeFem++  /Applications/FreeFem++.app/Contents/ff-4.12/share/FreeFEM/4.12/examples/tutorial/Laplace.edp 
   
if you have message the FreeFem++.app is corrupted then do in terminal
> sudo xattr -rc /Applications/FreeFem++.app
and retry. 
        
If you need to recomple  you can also install xcode and brew (cf. https://docs.brew.sh/Installation )    
and do 

	 brew install gcc 
	 