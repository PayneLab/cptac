**Things to know**



*   This is how we're storing our package version in a single location, accessible as needed: option 3 on [https://packaging.python.org/guides/single-sourcing-package-version/#single-sourcing-the-version](https://packaging.python.org/guides/single-sourcing-package-version/#single-sourcing-the-version)
*   Streaming downloads: Currently, the largest data file we have is around 60 MB, and most are significantly smaller. With these sizes, we have no problem downloading data files all at once. If, for some reason, we need to later stream the downloads, you can refer to the following web pages (as of July 2019) to edit the download_file function within "cptac/cptac/file_download.py":
    *   [https://stackoverflow.com/questions/16694907/download-large-file-in-python-with-requests](https://stackoverflow.com/questions/16694907/download-large-file-in-python-with-requests)
    *   [http://masnun.com/2016/09/18/python-using-the-requests-module-to-download-large-files-efficiently.html](http://masnun.com/2016/09/18/python-using-the-requests-module-to-download-large-files-efficiently.html)
    *   [https://www.geeksforgeeks.org/downloading-files-web-using-python/](https://www.geeksforgeeks.org/downloading-files-web-using-python/)
