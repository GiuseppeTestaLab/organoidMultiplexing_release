
import os
import argparse
import yaml
import markdown

INDEX_TEMPLATE_HEADER = """
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="X-UA-Compatible" content="ie=edge">
    <title>{1}</title>
    <link rel="stylesheet" href="./style.css">
    <link rel="icon" href="./favicon.ico" type="image/x-icon">
  </head>
  <body>
    <main>
        <h1>{0}</h1>
        {2}
    <p>
"""

INDEX_TEMPLATE_LINKLINE = "<h2><a href={1}>{0}</a></h2>{2}"

INDEX_TEMPLATE_FINAL = """
      {0}
    </p>
    </main>
    <script src="index.js"></script>
  </body>
</html>
"""

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("configyaml")
    args = parser.parse_args()
    
    with open(args.configyaml) as f:
      config_dict = yaml.load(f, Loader=yaml.SafeLoader)

    html_subtitle = markdown.markdown(config_dict["subtitle"])
    print(INDEX_TEMPLATE_HEADER.format(config_dict["title"], config_dict["short_title"], html_subtitle))
    for k, v in config_dict["pages"].items():
       html_description = markdown.markdown(v["description"])
       print(INDEX_TEMPLATE_LINKLINE.format(k, v["notebook"].replace("ipynb", "html"), html_description))
    html_final = markdown.markdown(config_dict["final"])
    print(INDEX_TEMPLATE_FINAL.format(html_final))

if __name__ == '__main__':
    main()
