export default {
  title: "C. elegans TE Explorer",
  root: "src",
  base: "/Low-input-ribosome-profiling-of-Celegans-embryogenesis/",
  theme: "air",
  sidebar: true,
  toc: false,
  pages: [
    {name: "Explorer", path: "/"},
    {name: "About", path: "/about"}
  ],
  head: `<style>
    :root {
      --theme-background: #fafafa;
    }
    .card {
      background: white;
      border: 1px solid #eaeaea;
      box-shadow: 0 2px 6px rgba(0,0,0,0.04);
      border-radius: 8px;
    }
    .card p {
      max-width: 100%;
    }
    #observablehq-main p,
    #observablehq-main h1,
    #observablehq-main h2,
    #observablehq-main h3,
    #observablehq-main table,
    #observablehq-main blockquote,
    #observablehq-main ul,
    #observablehq-main ol,
    #observablehq-main hr {
      max-width: 100%;
    }
  </style>`
};
