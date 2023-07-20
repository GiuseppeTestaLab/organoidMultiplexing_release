
import os
import argparse
import yaml

INDEX_TEMPLATE_HEADER = """
# {0}

**{1}**

{2}

An html version of the notebooks is accessible [here](https://{3}.github.io/{4}/).

"""

INDEX_TEMPLATE_LINKLINE = """

## {0}

Links: [jupyter notebook]({1}) and [html file]({2}).

{3}
"""

INDEX_TEMPLATE_LINKLINE_ONLY_HTML = """

## {0}

Links: [html file]({1}).

{2}
"""

INDEX_TEMPLATE_FINAL = """
{0}

---
*Note: this README file has been generated automatically.* <br>
*Please do not modify it directly but instead work on [this config file](resources/config.yaml).*

"""

EXCLUDED = ['index.html']

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("configyaml")
    args = parser.parse_args()
    
    with open(args.configyaml) as f:
      config_dict = yaml.load(f, Loader=yaml.SafeLoader)

    repo_full_name = os.environ['GITHUB_REPOSITORY']
    organization_name =  repo_full_name.split("/")[0]
    repo_name = repo_full_name.split("/")[1]
    
    print(INDEX_TEMPLATE_HEADER.format(config_dict["short_title"], config_dict["title"], config_dict["subtitle"],
                                       organization_name, repo_name))
    for k, v in config_dict["pages"].items():
      if os.path.splitext(v["notebook"])[1] == ".ipynb":
        ghpages_link = "https://{0}.github.io/{1}/{2}".format(organization_name, repo_name, v["notebook"].replace("ipynb", "html"))
        print(INDEX_TEMPLATE_LINKLINE.format(k, v["notebook"], ghpages_link, v["description"]))
      elif os.path.splitext(v["notebook"])[1] == ".html":
        ghpages_link = "https://{0}.github.io/{1}/{2}".format(organization_name, repo_name, v["notebook"])
        print(INDEX_TEMPLATE_LINKLINE_ONLY_HTML.format(k, ghpages_link, v["description"]))        
    print(INDEX_TEMPLATE_FINAL.format(config_dict["final"]))

if __name__ == '__main__':
    main()
